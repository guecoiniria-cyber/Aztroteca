#!/usr/bin/env python3
"""
Fichas científicas DAE (Dynamic Analog Engine) para las 10 zonas del ranking QGIS.

Entradas
  --resumen  DAE_resumen_analogía_preliminar_10_zonas.csv  (valores QGIS + puntuaciones)
  --matriz   DAE_matriz_62_variables_10_zonas.csv          (diccionario por zona)
  --qgis     carpeta qgis/ del proyecto (zonas_stats.gpkg, ranking_zonas.csv, rásters *_mx.tif)
  --salida   carpeta de salida (index.html + fichas/<fid>.html)

Lo que agrega a los datos del equipo (todo calculado, nada inventado):
  · estadísticas zonales por píxel (n, mediana, P5, P95) de precipitación, NDVI, amplitud y polvo
  · posición de cada zona frente a las 921 zonas del ranking nacional
  · DAE provisional de las 921 zonas → posición DAE nacional
  · sensibilidad del DAE a ±10 % en los pesos de precipitación y NDVI
  · confianza según el protocolo (proxy = 0.7; valor regional = 0.5)
  · extensión del polígono y latitud (variables de contexto/límite)
  · criterio de latitud del spaceport (0–28°)
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import math
import re
import sqlite3
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "fichas"))
import ficha_cientifica as fc  # noqa: E402  (reutiliza CSS, JS y lector GPKG)
sys.path.insert(0, str(Path(__file__).resolve().parent))
from lugares import LUGARES  # noqa: E402

# Capas oficiales mexicanas sugeridas (informe CONABIO) para cada variable pendiente
CAPAS_MX = {
    "Litología dominante": "CONABIO hidgeo4mgw (Descripcio contiene «basálticas y andesíticas») como filtro grueso; litología fina: SGM GeoInfoMex 1:250 000",
    "Mineralogía de alteración acuosa": "INEGI Edafología Serie II eda250s2gw, GRUPO1: Gypsisol→sulfatos, Calcisol→carbonatos, Vertisol→esmectitas, Solonchak→sales; EMIT para mineralogía directa",
    "Textura superficial": "Edafología INIFAP-CONABIO eda251mgw (textura y fase física)",
    "Formas y procesos eólicos": "Provincias fisiográficas rfisio4mgw («Desiertos»); Arenosol en eda250s2gw; topoformas INEGI",
    "Tipo de recurso hídrico": "eda250s2gw (Solonchak, Solonetz, fases salina/sódica); RHP rhpri4mgw",
    "Extremófilos": "Contexto: RTP rtp1mgw, Ramsar ramsar15gw y registros del SNIB",
    "Cobertura de rocas ≥0.5 m": "Proxy regional: Leptosol y fases líticas en eda250s2gw",
    "Pendiente local": "CEM INEGI 15 m (cálculo); DEM CONABIO demgw 1 km solo como contexto",
    "Pendiente transitable P95": "CEM INEGI 15 m (cálculo); DEM CONABIO demgw 1 km solo como contexto",
    "Temperatura máxima diurna de superficie": "Validación: Bioclimas bio05 (temperatura del aire, 1 km)",
    "Temperatura mínima nocturna de superficie": "Validación: Bioclimas bio06 (temperatura del aire, 1 km)",
    "Precipitación media": "Validación: Bioclimas bio12 ÷ 365 (CONABIO, CC BY 4.0, 1 km)",
    "Índice de aridez (P/ETP)": "Bioclimas bio12 con evapotranspiración evapr4mgw (CONABIO)",
    "Clima Köppen-García": "CONABIO clima1mgw, campo Clima_tipo (BW* árido, BS* semiárido)",
    "Uso de suelo y vegetación": "INEGI Serie VII usv250s7gw vía CONABIO, campo DESCRIPCION",
    "Área natural protegida": "CONANP anpenero2024gw vía CONABIO (unir por localización)",
    "Sitio Ramsar o región prioritaria": "ramsar15gw, rtp1mgw y rhpri4mgw (CONABIO)",
}


def clase_nac(n):
    return ("Análogo alto" if n >= 9 else "Análogo medio-alto" if n >= 7 else
            "Análogo moderado" if n >= 5 else "Análogo bajo")


def conformidad(f):
    """Rangos de la hoja del equipo «Análogos Tierra-Marte» (escala 1–5) con las variables que ya tienen dato."""
    la = abs(f["lat"]); pa = f["P"] * 365
    c = [
        ("Aridez / índice de vegetación (NDVI)", f"{f['Nd']:.3f}", "−0.1 a 0.15 → 5", 5 if -0.1 <= f["Nd"] <= 0.15 else 0),
        ("Tasa de precipitación real (GPM)", f"{pa:.0f} mm/año", "0 a 50 mm/año → 5", 5 if pa <= 50 else 0),
        ("Amplitud térmica diurna", f"{f['T']:.1f} °C", "25 a 250 °C → 5", 5 if 25 <= f["T"] <= 250 else 0),
        ("Latitud ecuatorial", f"{la:.2f}°", "0–28° → 5 · 29–50° → 4", 5 if la <= 28 else 4 if la <= 50 else 0),
        ("Extensión del terreno", f"{f['area_km2'] * 1e6:,.0f} m²", "≥ 747 000 m² → 5", 5 if f["area_km2"] * 1e6 >= 747000 else 0),
    ]
    return c, 100 * sum(x[3] for x in c) / (5 * len(c))


e, fmt = fc.e, fc.fmt

# ----------------------------------------------------------------------------------
# Parámetros del protocolo DAE (Diccionario de variables Tierra–Marte)
# ----------------------------------------------------------------------------------
W_P, W_N, W_T = 0.30, 0.30, 0.40          # 40 % provisional a la amplitud (única variable extra con dato)
RF_AMPLITUD, RF_RESTO = 6.0, 541.0        # cobertura = 60 % + 40 % · 6/541
Q_PROXY, Q_REGIONAL = 0.7, 0.5            # calidad del dato
MIN_PX = 4                                # píxeles mínimos para no ser "valor regional"
X0, Y0, D = -118.3745265, 32.7165358, 0.0083333   # georreferencia de los rásters (EPSG:4326, ~1 km)


def s_precip(p):
    return max(0.0, 1 - p / 2)


def s_ndvi(n):
    return 1.0 if n <= 0.1 else max(0.0, 1 - (n - 0.1) / 0.3)


def s_amp(t):
    if 40 <= t <= 100:
        return 1.0
    d = 40 - t if t < 40 else t - 100
    return max(0.0, 1 - d / 60)


def dae(p, n, t, wp=W_P, wn=W_N):
    wt = 1 - wp - wn
    return 100 * (wp * s_precip(p) + wn * s_ndvi(n) + wt * s_amp(t))


def estado(s):
    if s is None:
        return "sin dato"
    return "dentro" if s >= 0.999 else "cercano" if s >= 0.5 else "fuera"


def clase(i):
    return ("Análogo alto" if i >= 80 else "Análogo moderado" if i >= 60 else
            "Análogo parcial" if i >= 40 else "Análogo débil")


# ----------------------------------------------------------------------------------
# Lectura
# ----------------------------------------------------------------------------------
def leer_csv(p: Path):
    return list(csv.DictReader(io.StringIO(fc.leer_texto(p))))


def leer_matriz(p: Path):
    """Agrupa la matriz por FID y corrige filas con Variable/Rol intercambiados."""
    por_fid = {}
    for r in leer_csv(p):
        var, rol = r["Variable"].strip(), r["Rol"].strip()
        if var in ("Directa", "Proxy", "Categórica", "No puntúa") and rol not in ("Directa", "Proxy", "Categórica", "No puntúa"):
            var, rol = rol, var
        r["Variable"], r["Rol"] = var, rol
        por_fid.setdefault(int(r["FID"]), []).append(r)
    return por_fid


def leer_rasters(qdir: Path):
    import tifffile
    out = {}
    for clave, nombre in (("precip", "precip_mx.tif"), ("ndvi", "ndvi_limpio_mx.tif"),
                          ("rango", "rango_termico_mx.tif"), ("polvo", "polvo_mx.tif")):
        with tifffile.TiffFile(qdir / nombre) as t:
            a = t.pages[0].asarray().astype("float64")
        a[(a < -9000) | (a > 1e10)] = np.nan
        out[clave] = a
    return out


def zonal(a, polys):
    from matplotlib.path import Path as MP
    pts = [p for poly in polys for r in poly for p in r]
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    c0 = max(0, int((min(xs) - X0) / D)); c1 = min(a.shape[1], int((max(xs) - X0) / D) + 2)
    r0 = max(0, int((Y0 - max(ys)) / D)); r1 = min(a.shape[0], int((Y0 - min(ys)) / D) + 2)
    cc, rr = np.meshgrid(np.arange(c0, c1), np.arange(r0, r1))
    P = np.c_[(X0 + (cc + 0.5) * D).ravel(), (Y0 - (rr + 0.5) * D).ravel()]
    m = np.zeros(len(P), bool)
    for poly in polys:
        ins = MP(poly[0]).contains_points(P)
        for h in poly[1:]:
            ins &= ~MP(h).contains_points(P)
        m |= ins
    v = a[r0:r1, c0:c1].ravel()[m]
    v = v[~np.isnan(v)]
    if not len(v):
        return None
    return {"n": int(len(v)), "media": float(v.mean()), "mediana": float(np.median(v)),
            "p5": float(np.percentile(v, 5)), "p95": float(np.percentile(v, 95)),
            "min": float(v.min()), "max": float(v.max())}


def simplificar(polys, paso=0.02):
    out = []
    for poly in polys:
        anillos = []
        for r in poly:
            q, last = [], None
            for x, y in r:
                k = (round(x / paso) * paso, round(y / paso) * paso)
                if k != last:
                    q.append(k); last = k
            if len(q) >= 3:
                anillos.append(q)
        if anillos:
            out.append(anillos)
    return out


# ----------------------------------------------------------------------------------
# Cálculo por zona
# ----------------------------------------------------------------------------------
def pct_menor(valores, v):
    """% de zonas con valor estrictamente menor que v."""
    return 100 * float(np.mean(np.asarray(valores) < v))


def preparar(args):
    resumen = leer_csv(args.resumen)
    matriz = leer_matriz(args.matriz)
    nac = leer_csv(args.qgis / "ranking_zonas.csv")
    for z in nac:
        for k in ("ind_mean", "ndvi_mean", "lluvia_mean", "rango_mean", "polvo_mean"):
            z[k] = float(z[k])
        z["dae"] = dae(z["lluvia_mean"], z["ndvi_mean"], z["rango_mean"])
    orden = sorted(nac, key=lambda z: -z["dae"])
    rank_dae = {int(z["fid"]): i + 1 for i, z in enumerate(orden)}
    N = len(nac)
    arr = {k: np.array([z[k] for z in nac]) for k in ("ind_mean", "ndvi_mean", "lluvia_mean", "rango_mean", "polvo_mean", "dae")}

    con = sqlite3.connect(str(args.qgis / "zonas_stats.gpkg"))
    geoms = {r[0]: fc._gpkg_geom(r[1]) for r in con.execute("SELECT fid, geom FROM zonas_stats")}
    rasters = leer_rasters(args.qgis)

    nombres = {int(r["fid"]): r["name"] for r in resumen}
    fichas = []
    for r in resumen:
        fid = int(r["fid"])
        P, Nd, T, Pv = float(r["precip_raw"]), float(r["ndvi"]), float(r["thermal_range_raw"]), float(r["dust_raw"])
        sP, sN, sT = s_precip(P), s_ndvi(Nd), s_amp(T)
        d = dae(P, Nd, T)
        geo = fc.geometria_resumen(geoms[fid])
        zs = {k: zonal(a, geoms[fid]) for k, a in rasters.items()}
        regional = any(v and v["n"] < MIN_PX for v in zs.values())
        q = Q_REGIONAL if regional else Q_PROXY
        sens = [dae(P, Nd, T, wp, wn) for wp in (0.27, 0.33) for wn in (0.27, 0.33)]
        cob = 100 * (W_P + W_N + W_T * RF_AMPLITUD / RF_RESTO)
        lat = float(r["lat_centroid"]); lon = float(r["lon_centroid"])
        fichas.append({
            "fid": fid, "nombre": r["name"], "descripcion": r["description"],
            "rank_qgis": int(r["ranking_qgis"]), "rank_dae10": int(r["dae_rank_prelim"]),
            "rank_dae_nac": rank_dae[fid], "N": N,
            "lat": lat, "lon": lon, "area_km2": float(r["area_km2"]), "geo": geo,
            "qgis_index": float(r["qgis_index"]),
            "qgis_pct": pct_menor(arr["ind_mean"], float(r["qgis_index"])),
            "P": P, "Nd": Nd, "T": T, "Pv": Pv, "sP": sP, "sN": sN, "sT": sT,
            "dae": d, "dae_csv": float(r["DAE_prelim_pct"]),
            "aportes": {"Precipitación": 100 * W_P * sP, "NDVI": 100 * W_N * sN, "Amplitud térmica": 100 * W_T * sT},
            "cobertura": cob, "confianza": 100 * q, "regional": regional,
            "sens_min": min(sens), "sens_max": max(sens),
            "nac": {
                "precip": {"pct_mas_humedas": 100 - pct_menor(arr["lluvia_mean"], P) - 100 * float(np.mean(arr["lluvia_mean"] == P)),
                           "min": arr["lluvia_mean"].min(), "max": arr["lluvia_mean"].max(), "med": float(np.median(arr["lluvia_mean"]))},
                "ndvi": {"pct_mas_verdes": 100 * float(np.mean(arr["ndvi_mean"] > Nd)),
                         "min": arr["ndvi_mean"].min(), "max": arr["ndvi_mean"].max(), "med": float(np.median(arr["ndvi_mean"]))},
                "rango": {"pct_menor": pct_menor(arr["rango_mean"], T),
                          "min": arr["rango_mean"].min(), "max": arr["rango_mean"].max(), "med": float(np.median(arr["rango_mean"]))},
                "polvo": {"pct_menor": pct_menor(arr["polvo_mean"], Pv),
                          "min": arr["polvo_mean"].min(), "max": arr["polvo_mean"].max(), "med": float(np.median(arr["polvo_mean"]))},
                "dae": {"med": float(np.median(arr["dae"])), "max": float(arr["dae"].max())},
            },
            "zonal": zs, "matriz": matriz.get(fid, []), "nota": r["Nota_calculo"],
            "polys": geoms[fid],
        })
    for f in fichas:
        f["P_qgis"] = 100 * float(np.mean(arr["ind_mean"] < f["qgis_index"]))
        f["P_dae"] = 100 * float(np.mean(arr["dae"] < f["dae"]))
        f["nota"] = 1 + 9 * (f["P_qgis"] + f["P_dae"]) / 200
        f["nota_qgis"] = 1 + 9 * f["P_qgis"] / 100
        f["nota_dae"] = 1 + 9 * f["P_dae"] / 100
        f["clase_nac"] = clase_nac(f["nota"])
        f["conf_items"], f["conf"] = conformidad(f)
        f["lugares"] = [{"nombre": n, "lat": la, "lon": lo} for n, la, lo in LUGARES.get(f["fid"], [])]
    for k, f in enumerate(sorted(fichas, key=lambda f: (-round(f["nota"], 2), -f["dae"])), 1):
        f["rank_nota"] = k
    # zonas fuera del top-10 QGIS con DAE alto
    top_ids = {f["fid"] for f in fichas}
    otras = [{"fid": int(z["fid"]), "rank_qgis": int(z["ranking"]), "dae": z["dae"], "rank_dae": rank_dae[int(z["fid"])],
              "ndvi": z["ndvi_mean"], "precip": z["lluvia_mean"], "rango": z["rango_mean"]}
             for z in orden[:15] if int(z["fid"]) not in top_ids][:8]
    mapa = {"todas": {fid: simplificar(g) for fid, g in geoms.items()}, "top": top_ids,
            "puntos": {f["fid"]: (f["lon"], f["lat"]) for f in fichas},
            "lugares": {f["fid"]: f["lugares"] for f in fichas}}
    return fichas, otras, mapa, N


# ----------------------------------------------------------------------------------
# SVG: mapa de contexto y posición nacional
# ----------------------------------------------------------------------------------
def svg_mapa(mapa, resaltar=None, w=520, h=420, etiquetas=None, puntos=None):
    todas = mapa["todas"]
    xs = [x for g in todas.values() for p in g for r in p for x, _ in r]
    ys = [y for g in todas.values() for p in g for r in p for _, y in r]
    x0, x1, y0, y1 = min(xs) - .3, max(xs) + .3, min(ys) - .3, max(ys) + .3
    k = math.cos(math.radians((y0 + y1) / 2))
    sx = w / ((x1 - x0) * k); sy = h / (y1 - y0); s = min(sx, sy)
    ox = (w - (x1 - x0) * k * s) / 2; oy = (h - (y1 - y0) * s) / 2

    def path(g):
        return "".join("M" + "L".join(f"{ox + (x - x0) * k * s:.1f},{oy + (y1 - y) * s:.1f}" for x, y in r) + "Z"
                       for p in g for r in p)
    out = [f"<svg viewBox='0 0 {w} {h}' class='map' role='img' aria-label='Mapa de zonas candidatas'>"]
    out.append("<path class='z-all' d='" + "".join(path(g) for fid, g in todas.items() if fid not in mapa["top"]) + "'/>")
    out.append("<path class='z-top' d='" + "".join(path(todas[f]) for f in mapa["top"] if f != resaltar) + "'/>")
    if resaltar is not None:
        g = todas[resaltar]
        out.append(f"<path class='z-hi' d='{path(g)}'/>")
        lo_, la_ = puntos[resaltar]
        X = ox + (lo_ - x0) * k * s; Y = oy + (y1 - la_) * s
        out.append(f"<circle class='z-ring' cx='{X:.1f}' cy='{Y:.1f}' r='14'/>")
        for L in mapa.get("lugares", {}).get(resaltar, []):
            Xl = ox + (L["lon"] - x0) * k * s; Yl = oy + (y1 - L["lat"]) * s
            out.append(f"<circle class='z-pl' cx='{Xl:.1f}' cy='{Yl:.1f}' r='3'><title>{e(L['nombre'])}</title></circle>")
    if etiquetas:
        for fid, txt in etiquetas.items():
            lo_, la_ = puntos[fid]
            X = ox + (lo_ - x0) * k * s; Y = oy + (y1 - la_) * s
            out.append(f"<g class='lbl'><circle cx='{X:.1f}' cy='{Y:.1f}' r='9'/><text x='{X:.1f}' y='{Y + 3.5:.1f}'>{e(txt)}</text></g>")
    out.append("</svg>")
    return "".join(out)


def strip(valor, lo, hi, med, ventana=None, invertir=False, fmtv=lambda v: fmt(v), etiqueta=""):
    """Barra horizontal con el rango nacional (921 zonas), la mediana, el sitio y la ventana marciana."""
    span = (hi - lo) or 1

    def pos(v):
        return max(0, min(100, (v - lo) / span * 100))
    H = ["<div class='strip'>"]
    if ventana:
        a, b = pos(max(lo, ventana[0])), pos(min(hi, ventana[1]))
        if ventana[1] >= lo and ventana[0] <= hi and b > a:
            H.append(f"<span class='win' style='left:{a:.1f}%;width:{max(1.2, b - a):.1f}%' title='Ventana marciana'></span>")
    H.append(f"<span class='med' style='left:{pos(med):.1f}%' title='Mediana nacional'></span>")
    H.append(f"<span class='me' style='left:{pos(valor):.1f}%'><b>{fmtv(valor)}</b></span>")
    H.append(f"</div><div class='strip-ax'><span>{fmtv(lo)}</span><span>{e(etiqueta)}</span><span>{fmtv(hi)}</span></div>")
    return "".join(H)


# ----------------------------------------------------------------------------------
# HTML de la ficha
# ----------------------------------------------------------------------------------
EXTRA_CSS = r"""
.hero .desc{font-size:17px;color:var(--ink2);max-width:640px;margin:-6px 0 20px;line-height:1.55}
.prov{font:600 10.5px var(--mono);letter-spacing:.14em;color:var(--mid);text-transform:uppercase;margin-top:10px}
.idx{display:grid;grid-template-columns:repeat(3,1fr);gap:14px}
.ic{padding:24px;position:relative}.ic .l{font:600 11px var(--mono);letter-spacing:.16em;text-transform:uppercase;color:var(--mut)}
.ic .v{font:700 40px/1.1 var(--disp);letter-spacing:-.02em;margin:10px 0 6px}.ic .v small{font-size:16px;color:var(--mut);font-weight:500;margin-left:3px}
.ic p{font-size:13px;color:var(--ink2);margin:10px 0 0;line-height:1.55}.ic .f{font:500 12px var(--mono);color:var(--mut);margin-top:10px;padding-top:10px;border-top:1px solid var(--line)}
.ic.wide{grid-column:span 2}
.stack{display:flex;height:34px;border-radius:10px;overflow:hidden;background:rgba(255,255,255,.05);margin:16px 0 10px}
.stack i{display:flex;align-items:center;justify-content:center;font:600 12px var(--mono);color:#fff;width:0;transition:width 1.4s cubic-bezier(.2,.8,.2,1);white-space:nowrap;overflow:hidden}
.stack i:nth-child(1){background:#c1440e}.stack i:nth-child(2){background:#e0782f}.stack i:nth-child(3){background:#8a5a3c}
.leg{display:flex;flex-wrap:wrap;gap:8px 18px;font-size:12.5px;color:var(--ink2)}.leg span{display:flex;align-items:center;gap:7px}
.leg i{width:10px;height:10px;border-radius:3px}
.range{position:relative;height:8px;background:rgba(255,255,255,.07);border-radius:9px;margin:18px 0 8px}
.range .r{position:absolute;top:0;bottom:0;background:rgba(242,177,52,.35);border-radius:9px}
.range .p{position:absolute;top:-5px;width:3px;height:18px;background:var(--ink);border-radius:2px}
.nat{display:grid;grid-template-columns:repeat(2,1fr);gap:14px}
.nc{padding:22px 24px}.nc .t{display:flex;justify-content:space-between;gap:12px;align-items:baseline}
.nc .t b{font:600 16px var(--disp)}.nc .t span{font:500 12px var(--mono);color:var(--mut)}
.nc .big{font:700 30px/1.1 var(--disp);margin:10px 0 2px}.nc .big small{font-size:14px;color:var(--mut);font-weight:500;margin-left:4px}
.nc p{font-size:13px;color:var(--ink2);margin:6px 0 0}
.strip{position:relative;height:12px;margin:34px 0 6px;background:linear-gradient(90deg,rgba(255,255,255,.04),rgba(255,255,255,.1));border-radius:9px}
.strip .win{position:absolute;top:0;bottom:0;background:rgba(62,207,142,.28);border:1px solid rgba(62,207,142,.6);border-radius:9px}
.strip .med{position:absolute;top:-4px;bottom:-4px;width:2px;background:var(--mut);transform:translateX(-1px)}
.strip .me{position:absolute;top:-6px;width:14px;height:24px;transform:translateX(-7px);background:var(--burn);border-radius:5px;box-shadow:0 0 12px var(--glow)}
.strip .me b{position:absolute;bottom:30px;left:50%;transform:translateX(-50%);font:600 12px var(--mono);white-space:nowrap;color:var(--ink)}
.strip-ax{display:flex;justify-content:space-between;font:500 11px var(--mono);color:var(--mut2)}
.slegend{display:flex;flex-wrap:wrap;gap:6px 16px;font-size:12px;color:var(--mut);margin-top:14px}.slegend span{display:flex;gap:6px;align-items:center}
.slegend i{display:inline-block;width:12px;height:8px;border-radius:3px}
.mapcard{padding:18px}.map{width:100%;height:auto;display:block}
.z-all{fill:rgba(255,255,255,.10);stroke:rgba(255,255,255,.18);stroke-width:.3}
.z-top{fill:rgba(224,120,47,.35);stroke:#e0782f;stroke-width:.5}
.z-hi{fill:#c1440e;stroke:#ffd9c4;stroke-width:.9;filter:drop-shadow(0 0 6px rgba(193,68,14,.9))}
.z-ring{fill:none;stroke:#ffd9c4;stroke-width:1.2;opacity:.8;animation:ring 2.4s infinite}
.z-dot{fill:#ffd9c4}@keyframes ring{0%{r:6;opacity:.9}100%{r:22;opacity:0}}
.lbl circle{fill:#c1440e;stroke:#ffd9c4;stroke-width:1}.lbl text{font:700 9px var(--mono);fill:#fff;text-anchor:middle}
.mapcap{font-size:12px;color:var(--mut);margin-top:10px;display:flex;gap:14px;flex-wrap:wrap}.mapcap i{display:inline-block;width:10px;height:10px;border-radius:3px;margin-right:6px;vertical-align:-1px}
.zt{width:100%;border-collapse:collapse;font-size:13.5px}
.zt th{font:600 10.5px var(--mono);letter-spacing:.12em;text-transform:uppercase;color:var(--mut2);text-align:left;padding:12px 16px;border-bottom:1px solid var(--line);background:var(--bg2)}
.zt td{padding:12px 16px;border-bottom:1px solid var(--line);font-family:var(--mono);font-size:13px}
.zt td:first-child{font-family:var(--body);font-weight:600}.zt tr:hover td{background:rgba(255,255,255,.025)}
.callout{padding:20px 24px;border-left:3px solid var(--mid);background:rgba(242,177,52,.06);border-radius:0 12px 12px 0;color:var(--ink2);font-size:14px;margin-top:16px}
.callout b{color:var(--ink)}
.rk{display:grid;grid-template-columns:44px 1.6fr 1fr 110px 90px 90px;gap:16px;align-items:center;padding:16px 22px;border-top:1px solid var(--line);text-decoration:none;transition:background .2s}
.rk:hover{background:rgba(255,255,255,.03)}.rk:first-child{border-top:0}
.rk .n{font:700 22px var(--disp);color:var(--burn)}.rk .t b{display:block;font:600 15.5px var(--disp)}.rk .t span{font-size:12.5px;color:var(--mut)}
.rk .v{font:700 20px var(--disp);text-align:right}.rk .m{font:500 12.5px var(--mono);color:var(--ink2);text-align:right}
.rk.hd{font:600 10.5px var(--mono);letter-spacing:.12em;text-transform:uppercase;color:var(--mut2);background:var(--bg2);padding:12px 22px}
.rk.hd span{text-align:right}
.z-pl{fill:#fff;stroke:#c1440e;stroke-width:1.2}
.kpis.k5{grid-template-columns:repeat(5,1fr)}
.scale{position:relative;margin:26px 0 8px}
.scale .bands{display:flex;height:14px;border-radius:9px;overflow:hidden}
.scale .bands i{display:block;height:100%}
.scale .mk{position:absolute;top:-8px;width:4px;height:30px;background:#fff;border-radius:2px;transform:translateX(-2px);box-shadow:0 0 12px rgba(255,255,255,.6)}
.scale .mk b{position:absolute;bottom:34px;left:50%;transform:translateX(-50%);font:700 14px var(--mono);white-space:nowrap}
.scale-ax{display:flex;justify-content:space-between;font:500 11px var(--mono);color:var(--mut2);margin-top:6px}
.scale-lg{display:flex;flex-wrap:wrap;gap:6px 16px;font-size:12px;color:var(--ink2);margin-top:12px}.scale-lg span{display:flex;gap:6px;align-items:center}.scale-lg i{width:10px;height:10px;border-radius:3px}
.pt{font:600 12px var(--mono)}.pt.ok{color:var(--ok)}.pt.bad{color:var(--bad)}.pt.mid{color:var(--mid)}
.scat{width:100%;height:auto;display:block}.scat text{font:500 11px var(--mono);fill:var(--mut)}
.scat .ax{stroke:var(--line2)}.scat .gd{stroke:rgba(255,255,255,.05)}.scat .pt-c{fill:#c1440e;stroke:#ffd9c4;stroke-width:1}
.scat .lbt{font:600 11.5px var(--body);fill:var(--ink)}.scat .band{fill:rgba(62,207,142,.06)}
@media (max-width:980px){.kpis.k5{grid-template-columns:repeat(2,1fr)}}.rk.hd span:nth-child(2){text-align:left}
@media (max-width:980px){.idx,.nat{grid-template-columns:1fr}.ic.wide{grid-column:auto}.rk{grid-template-columns:36px 1fr 80px}.rk .bar-c,.rk .m{display:none}.rk.hd{display:none}}
"""


def pct(v, nd=2):
    return f"{v:.{nd}f}".rstrip("0").rstrip(".") if nd else f"{v:.0f}"


def fila_var(r, valor, s, nota, nombre_tip):
    est = estado(s)
    return (f"<div class='vr' data-estado='{e(est)}'><div class='vn'><span>{e(r['Variable'])}</span>{fc._tip(nombre_tip)}</div>"
            f"<div class='c-val'><div class='val'>{e(valor)}</div>" + (f"<div class='note'>{e(nota)}</div>" if nota else "") + "</div>"
            f"<div class='c-rng rng'>{e(r['Referencia_Marte'])}</div><div>{fc._pill(est)}</div>"
            f"<div class='sco'>{(f'{s*100:.1f}%' if s is not None else '—')}</div>"
            f"<div class='rf'>{e(r['Peso_RxF'] or ('30 % fijo' if r['Variable'] in ('Precipitación media', 'NDVI') else '—'))}</div>"
            "<div><button type='button' class='srcb' aria-expanded='false' title='Método y fuente'>SRC</button></div>"
            f"<div class='src'><b>Método / fuente México:</b> {e(r['Metodo_o_fuente'])}<br>"
            + (f"<b>Capa oficial mexicana:</b> {e(CAPAS_MX[r['Variable']])}<br>" if r['Variable'] in CAPAS_MX else "") +
            f"<b>Rol:</b> {e(r['Rol'])} · <b>Estado en la matriz:</b> {e(r['Unidad_o_estado'])}</div></div>")


POR_QUE = {
    "Precipitación media": "Marte no tiene agua líquida estable en superficie; cuanto más cerca de 0 mm/día, más se parece el sitio a la hiperaridez marciana. Pesa 30 % fijo del índice.",
    "NDVI": "Marte no tiene vegetación; NDVI ≤ 0.1 indica roca, arena o costra desnuda. Pesa 30 % fijo del índice.",
    "Amplitud térmica diaria": "Las oscilaciones térmicas fuertes fatigan sellos, rocas y estructuras; en Marte superan 40 °C diarios. Es un proxy MODIS (fidelidad 2).",
}


def escala_html(nota):
    pos = (nota - 1) / 9 * 100
    return ("<div class='scale'><div class='bands'>"
            "<i style='width:44.44%;background:rgba(240,87,63,.45)'></i><i style='width:22.22%;background:rgba(242,177,52,.45)'></i>"
            "<i style='width:22.22%;background:rgba(224,120,47,.55)'></i><i style='width:11.12%;background:rgba(62,207,142,.65)'></i></div>"
            f"<span class='mk' style='left:{pos:.1f}%'><b>{nota:.2f}</b></span></div>"
            "<div class='scale-ax'><span>1</span><span>5</span><span>7</span><span>9</span><span>10</span></div>"
            "<div class='scale-lg'><span><i style='background:rgba(240,87,63,.6)'></i>Bajo &lt; 5</span><span><i style='background:rgba(242,177,52,.6)'></i>Moderado 5–7</span>"
            "<span><i style='background:rgba(224,120,47,.7)'></i>Medio-alto 7–9</span><span><i style='background:rgba(62,207,142,.8)'></i>Alto ≥ 9</span></div>")


def ficha_html(f, mapa, nombres):
    i = f["dae"]; g = f["geo"]; zs = f["zonal"]; n = f["nac"]
    mat = f["matriz"]
    pend_cnt = {"dentro": 0, "cercano": 0, "fuera": 0, "sin dato": 0}
    vals = {"Precipitación media": (f"{f['P']:.3f} mm/día", f["sP"], "Unidad asumida mm/día (GPM IMERG); validar conversión y periodo."),
            "NDVI": (f"{f['Nd']:.4f}", f["sN"], "MOD13A3; validar escala ×0.0001 y periodo."),
            "Amplitud térmica diaria": (f"{f['T']:.2f} Δ°C", f["sT"], f"Por debajo del mínimo marciano (40 °C) por {40 - f['T']:.1f} °C. Proxy MODIS día − noche.")}
    puntuables = [r for r in mat if r["Rol"] in ("Directa", "Proxy", "Categórica")]
    for r in puntuables:
        s = vals.get(r["Variable"], (None, None, None))[1]
        pend_cnt[estado(s)] += 1
    clas = clase(i)

    H = ["<!DOCTYPE html><html lang='es'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>",
         f"<title>{e(f['nombre'])} · Ficha DAE de análogo marciano</title>",
         f"<style>{fc.FICHA_CSS}{EXTRA_CSS}</style></head><body>"]
    nav = [("resumen", "Resumen"), ("indices", "Índices"), ("nacional", "Contexto nacional"), ("balance", "Balance"),
           ("variables", "Variables"), ("limites", "Límites"), ("spaceport", "Spaceport"), ("metodo", "Método")]
    H.append(f"<header class='top'><div class='wrap'><a class='brand' href='../index.html' style='text-decoration:none'><i></i>DAE · FID {f['fid']}</a><nav class='nav'>"
             + "".join(f"<a href='#{a}'>{t}</a>" for a, t in nav) + "</nav></div></header>")

    # HERO
    H.append("<section class='hero' id='top'><div class='wrap'><div class='hgrid'><div>"
             f"<div class='kick'>Ficha científica · #{f['rank_nota']} de 10 en la escala DAE-MX · #{f['rank_qgis']} QGIS</div>"
             f"<h1>{e(f['nombre'])}</h1><p class='desc'>{e(f['descripcion'])}</p><div class='chips'>"
             f"<span class='chip'>FID {f['fid']}</span><span class='chip'>{f['lat']:.4f}°, {f['lon']:.4f}°</span>"
             f"<span class='chip'>Área {f['area_km2']:,.1f} km²</span>"
             + (f"<span class='chip'>Extensión {g['extension_km']:,.1f} km</span>" if g else "") +
             f"<span class='chip'>{zs['rango']['n'] if zs['rango'] else '—'} píxeles de ~1 km</span></div>"
             f"<div class='cls'><i></i>{e(f['clase_nac'])} en México · nota {f['nota']:.2f}/10</div></div>")
    gv = (f["nota"] - 1) / 9
    H.append("<div class='gauge'><svg viewBox='0 0 120 120'><defs><linearGradient id='gg' x1='0' y1='0' x2='1' y2='1'>"
             "<stop offset='0' stop-color='#e0782f'/><stop offset='1' stop-color='#c1440e'/></linearGradient></defs>"
             "<circle class='trk' cx='60' cy='60' r='54' fill='none' stroke-width='7'/>"
             f"<circle class='val' data-v='{gv:.4f}' cx='60' cy='60' r='54' fill='none' stroke-width='7'/></svg>"
             f"<div class='ctr'><div class='big' data-count='{f['nota']:.2f}' data-dec='2'>{f['nota']:.2f}</div>"
             "<div class='lab'>Nota de analogía · México</div><div class='sub'>escala 1–10 · 921 zonas</div></div></div></div>")
    kp = [("Nota nacional (1–10)", f"<span data-count='{f['nota']:.2f}' data-dec='2'>{f['nota']:.2f}</span>", f"{f['clase_nac']} · #{f['rank_nota']} de 10", gv),
          ("DAE frente a Marte", f"<span data-count='{i:.2f}' data-dec='2'>{i:.2f}</span><small>%</small>", f"{clas} · #{f['rank_dae_nac']} de {f['N']}", i / 100),
          ("Conformidad rangos equipo", f"<span data-count='{f['conf']:.0f}'>{f['conf']:.0f}</span><small>%</small>", "5 de 45 criterios con dato", f["conf"] / 100),
          ("Cobertura de datos", f"<span data-count='{f['cobertura']:.2f}' data-dec='2'>{f['cobertura']:.2f}</span><small>%</small>", "Peso metodológico con dato", f["cobertura"] / 100),
          ("Confianza", f"<span data-count='{f['confianza']:.0f}'>{f['confianza']:.0f}</span><small>%</small>", "3 proxies (calidad 0.7)" + (" · valor regional" if f["regional"] else ""), f["confianza"] / 100)]
    H.append("<div class='kpis k5'>")
    for t, v, d, b in kp:
        H.append(f"<div class='kpi'><div class='l'>{t}</div><div class='v'>{v}</div><div class='d'>{e(d)}</div>{fc._bar(b, '')}</div>")
    H.append("</div><div class='tally'>")
    for est in ("dentro", "cercano", "fuera", "sin dato"):
        H.append(f"<button type='button' class='tbtn' data-f='{est}'><span class='dot'></span><b>{pend_cnt[est]}</b>{fc.ESTADO_TXT[est]}</button>")
    H.append("</div></div></section><main class='wrap'>")

    # RESUMEN
    fuertes = [x for x, s in (("NDVI", f["sN"]), ("Precipitación", f["sP"])) if s >= 0.8]
    p1 = (f"{f['nombre']} tiene nota {f['nota']:.2f} de 10 en la escala nacional: {f['clase_nac'].lower()} dentro de México. "
          f"Frente a Marte su DAE provisional es {i:.2f} % ({clas.lower()}), con {f['cobertura']:.2f} % del peso metodológico respaldado por datos. "
          f"Por DAE es la #{f['rank_dae10']} de las 10 zonas del top QGIS (donde ocupa el #{f['rank_qgis']}) y la #{f['rank_dae_nac']} de las {f['N']} zonas candidatas del país.")
    veg = ("Su vegetación es casi nula" if f['Nd'] <= 0.13 else "Conserva algo de vegetación")
    rel = ('dentro de' if f['sN'] >= .999 else 'cerca de' if f['sN'] >= .8 else 'lejos de')
    p2 = (f"{veg} (NDVI {f['Nd']:.4f}; {rel} la ventana marciana ≤ 0.1) "
          f"y llueve {f['P']:.2f} mm/día: solo el {100 - n['precip']['pct_mas_humedas']:.0f} % de las zonas candidatas es más seco. ")
    p3 = (f"La principal diferencia con Marte es la amplitud térmica: {f['T']:.1f} °C frente al mínimo marciano de 40 °C. "
          f"Ese componente aporta {f['aportes']['Amplitud térmica']:.1f} de los 40 puntos posibles, y es un proxy MODIS que tiende a subestimar la amplitud real.")
    if f["Pv"] >= n["polvo"]["med"] * 2:
        p4 = f"Tiene más polvo atmosférico que el {n['polvo']['pct_menor']:.0f} % de las zonas, lo que la hace útil para ensayos de mitigación de polvo (variable aún sin puntuar)."
    else:
        p4 = f"Su carga de polvo atmosférico supera a la del {n['polvo']['pct_menor']:.0f} % de las zonas; la variable aún no puntúa hasta confirmar la unidad de MERRA-2."
    usos = ["Ensayos de operación en hiperaridez: equipos, paneles solares y gestión de agua"]
    if f["sN"] >= .9:
        usos.append("Pruebas de movilidad y navegación sobre superficie desnuda")
    if f["Pv"] >= n["polvo"]["med"] * 2:
        usos.append("Pruebas de mitigación de polvo en paneles, sellos y óptica")
    usos.append("Caracterización en campo de regolito, geología y mineralogía para completar el índice")
    H.append("<section class='blk rv' id='resumen'><div class='sh'><div><div class='n'>01 · Resumen</div><h2>¿Por qué es análogo de Marte?</h2></div></div>"
             "<div class='why'><div class='card story'>"
             f"<p>{e(p1)}</p><p>{e(p2)}</p><p>{e(p3)}</p><p>{e(p4)}</p></div>"
             "<div class='card mapcard'>" + svg_mapa(mapa, resaltar=f["fid"], puntos=mapa['puntos']) +
             "<div class='mapcap'><span><i style='background:#c1440e'></i>Esta zona</span><span><i style='background:rgba(224,120,47,.5)'></i>Top 10 QGIS</span>"
             f"<span><i style='background:rgba(255,255,255,.18)'></i>{f['N']} zonas candidatas</span></div></div></div>"
             "<div class='card usos' style='margin-top:16px'><h3>Usos recomendados del análogo (según los datos actuales)</h3>"
             + "".join(f"<div class='uso'><span>{k:02d}</span><div>{e(u)}</div></div>" for k, u in enumerate(usos, 1)) + "</div>"
             "<div class='sub2'><h3>Lugares de referencia</h3><span>documento de zonas del equipo</span></div><div class='card list'>"
             "<div class='lr hd'><span>Lugar</span><span>Coordenadas</span><span>Verificación con el polígono QGIS</span></div>"
             + "".join(f"<div class='lr'><span class='k'>{e(L['nombre'])}</span><span class='m'>{L['lat']:.3f}°, {L['lon']:.3f}°</span>"
                       f"<span class='s'><span class='pt ok'>✓ dentro del polígono</span></span></div>" for L in f["lugares"])
             + "</div></section>")

    # ÍNDICES
    ap = f["aportes"]
    H.append("<section class='blk rv' id='indices'><div class='sh'><div><div class='n'>02 · Índices</div><h2>Cómo se forma cada porcentaje</h2></div>"
             "<p>Tres resultados separados: analogía, cobertura y confianza. El índice QGIS se muestra aparte porque usa otros criterios.</p></div>"
             "<div class='idx'>"
             f"<div class='card ic wide'><div class='l'>Escala de analogía nacional DAE-MX (1–10)</div><div class='v'>{f['nota']:.2f}<small>/ 10 · {e(f['clase_nac'])}</small></div>"
             + escala_html(f["nota"]) +
             f"<p>Promedio del percentil de la zona en dos índices independientes, entre las {f['N']} zonas candidatas: índice QGIS (supera al {f['P_qgis']:.1f} %) y DAE (supera al {f['P_dae']:.1f} %). "
             f"Con solo QGIS la nota sería {f['nota_qgis']:.2f}; con solo DAE, {f['nota_dae']:.2f}.</p>"
             "<div class='f'>Nota = 1 + 9 × (P_QGIS + P_DAE) / 200. Alto ≥ 9 equivale al ~11 % superior del país. Mide analogía relativa dentro de México; la analogía frente a Marte es el DAE.</div></div>"
             f"<div class='card ic'><div class='l'>Conformidad con rangos del equipo</div><div class='v'>{f['conf']:.0f}<small>%</small></div>{fc._bar(f['conf'] / 100, '')}"
             "<p>Hoja «Rangos de Variables» (Tierra–Marte, escala 1–5). Solo 5 de sus 45 criterios tienen dato todavía.</p>"
             "<div class='f'>" + "<br>".join(f"{e(a)}: {e(b)} → <span class='pt {'ok' if v == 5 else 'mid' if v else 'bad'}'>{v}/5</span>" for a, b, c, v in f["conf_items"]) + "</div></div>"
             f"<div class='card ic wide'><div class='l'>DAE provisional · aporte por componente</div><div class='v'>{i:.2f}<small>de 100 puntos</small></div>"
             "<div class='stack'>" + "".join(f"<i data-w='{v:.2f}' title='{e(k)}: {v:.2f} puntos'>{v:.1f}</i>" for k, v in ap.items()) + "</div>"
             "<div class='leg'>"
             f"<span><i style='background:#c1440e'></i>Precipitación 0.30 × {f['sP']:.4f} = {ap['Precipitación']:.2f}</span>"
             f"<span><i style='background:#e0782f'></i>NDVI 0.30 × {f['sN']:.4f} = {ap['NDVI']:.2f}</span>"
             f"<span><i style='background:#8a5a3c'></i>Amplitud 0.40 × {f['sT']:.4f} = {ap['Amplitud térmica']:.2f}</span></div>"
             f"<div class='f'>DAE* = 100 (0.30 s_P + 0.30 s_N + 0.40 s_T). El 40 % va provisionalmente a la amplitud por ser la única variable extra con dato. Coincide con el CSV del equipo: {f['dae_csv']:.2f} %.</div></div>"
             f"<div class='card ic'><div class='l'>Sensibilidad de pesos ±10 %</div><div class='v'>{f['sens_min']:.1f}–{f['sens_max']:.1f}<small>%</small></div>"
             f"<div class='range'><span class='r' style='left:{f['sens_min']:.1f}%;width:{max(.8, f['sens_max'] - f['sens_min']):.1f}%'></span><span class='p' style='left:{i:.1f}%'></span></div>"
             "<p>Rango del DAE si los pesos de precipitación y NDVI varían entre 27 % y 33 %, con el resto para la amplitud (protocolo, paso 10).</p></div>"
             f"<div class='card ic'><div class='l'>Cobertura de datos</div><div class='v'>{f['cobertura']:.2f}<small>%</small></div>{fc._bar(f['cobertura'] / 100, '')}"
             "<p>Parte del peso metodológico que ya tiene dato. No es parecido con Marte.</p><div class='f'>60 % + 40 % × (6 / 541) — R×F de la amplitud sobre el total R×F restante.</div></div>"
             f"<div class='card ic'><div class='l'>Confianza</div><div class='v'>{f['confianza']:.0f}<small>%</small></div>{fc._bar(f['confianza'] / 100, '')}"
             f"<p>Las 3 variables con dato son proxies (calidad 0.7). Cada una tiene {min(v['n'] for v in zs.values() if v)}+ píxeles, así que no hay valores regionales.</p>"
             "<div class='f'>Campo o satélite directo = 1 · proxy = 0.7 · valor regional = 0.5</div></div>"
             f"<div class='card ic'><div class='l'>Índice QGIS (ind_mean)</div><div class='v'>{f['qgis_index']:.3f}</div>{fc._bar((f['qgis_index'] - 0.73) / (0.81 - 0.73), '')}"
             f"<p>Posición #{f['rank_qgis']} de {f['N']}: supera al {f['qgis_pct']:.1f} % de las zonas. Índice zonal de QGIS con 4 componentes (precipitación, NDVI, amplitud y polvo).</p>"
             "<div class='f'>No equivale al DAE: se construyó con otros criterios.</div></div>"
             "</div></section>")

    # CONTEXTO NACIONAL
    H.append("<section class='blk rv' id='nacional'><div class='sh'><div><div class='n'>03 · Contexto nacional</div><h2>Posición frente a las 921 zonas candidatas</h2></div>"
             "<p>Cada barra va del mínimo al máximo nacional. Marca naranja = esta zona; línea gris = mediana; verde = ventana marciana.</p></div><div class='nat'>")
    cards = [
        ("Precipitación media", f"{f['P']:.3f}", "mm/día", f"Más seca que el {n['precip']['pct_mas_humedas']:.0f} % de las zonas.",
         strip(f["P"], n["precip"]["min"], n["precip"]["max"], n["precip"]["med"], (0, 0.012), fmtv=lambda v: f"{v:.2f}", etiqueta="mm/día")),
        ("NDVI", f"{f['Nd']:.4f}", "", f"{n['ndvi']['pct_mas_verdes']:.0f} % de las zonas tiene más vegetación.",
         strip(f["Nd"], n["ndvi"]["min"], n["ndvi"]["max"], n["ndvi"]["med"], (0, 0.1), fmtv=lambda v: f"{v:.3f}", etiqueta="NDVI")),
        ("Amplitud térmica", f"{f['T']:.2f}", "Δ°C", f"Mayor que en el {n['rango']['pct_menor']:.0f} % de las zonas; Marte pide ≥ 40 °C (fuera de escala).",
         strip(f["T"], n["rango"]["min"], n["rango"]["max"], n["rango"]["med"], None, fmtv=lambda v: f"{v:.1f}", etiqueta="Δ°C")),
        ("Polvo atmosférico (sin puntuar)", f"{f['Pv'] * 1e9:.1f}", "µg/m³ si es DUSMASS", f"Más polvo que el {n['polvo']['pct_menor']:.0f} % de las zonas.",
         strip(f["Pv"] * 1e9, n["polvo"]["min"] * 1e9, n["polvo"]["max"] * 1e9, n["polvo"]["med"] * 1e9, None, fmtv=lambda v: f"{v:.1f}", etiqueta="µg/m³*")),
    ]
    for t, v, u, p, s in cards:
        H.append(f"<div class='card nc'><div class='t'><b>{e(t)}</b><span>{e(u)}</span></div><div class='big'>{v}<small>{e(u.split(' ')[0])}</small></div><p>{e(p)}</p>{s}</div>")
    H.append("</div><div class='slegend'><span><i style='background:var(--burn)'></i>Esta zona</span><span><i style='background:var(--mut);width:2px'></i>Mediana nacional</span>"
             "<span><i style='background:rgba(62,207,142,.4)'></i>Ventana marciana</span><span>* Polvo: valor MERRA-2 ×10⁹; unidad por confirmar.</span></div>")
    # variabilidad interna
    H.append("<div class='sub2'><h3>Variabilidad dentro del polígono</h3><span>estadística zonal por píxel</span></div><div class='card list' style='overflow-x:auto'>"
             "<table class='zt'><tr><th>Variable</th><th>Píxeles</th><th>Media</th><th>Mediana</th><th>P5</th><th>P95</th><th>Mín.</th><th>Máx.</th></tr>")
    for lab, k, mult, nd in (("Precipitación (mm/día)", "precip", 1, 3), ("NDVI", "ndvi", 1, 4), ("Amplitud térmica (Δ°C)", "rango", 1, 2), ("Polvo (×10⁻⁹)", "polvo", 1e9, 2)):
        z = zs[k]
        if z:
            H.append(f"<tr><td>{lab}</td><td>{z['n']}</td>" + "".join(f"<td>{z[c] * mult:.{nd}f}</td>" for c in ("media", "mediana", "p5", "p95", "min", "max")) + "</tr>")
    H.append("</table></div><p class='note' style='color:var(--mut);font-size:13px;margin-top:10px'>La media coincide con el valor zonal de QGIS. El protocolo DAE recomienda usar la mediana; aquí se conserva la media para no cambiar el ranking del equipo.</p></section>")

    # BALANCE
    comps = [("NDVI", f"{f['Nd']:.4f}", f["sN"]), ("Precipitación media", f"{f['P']:.3f} mm/día", f["sP"]), ("Amplitud térmica", f"{f['T']:.2f} Δ°C", f["sT"])]

    def lst(cond, vacio):
        it = [(a, b, s) for a, b, s in comps if cond(s)]
        if not it:
            return f"<div class='empty'>{vacio}</div>"
        return "<ul>" + "".join(f"<li><span>{e(a)} · {s * 100:.1f}%</span><b>{e(b)}</b></li>" for a, b, s in it) + "</ul>"
    H.append("<section class='blk rv' id='balance'><div class='sh'><div><div class='n'>04 · Balance</div><h2>Fortalezas y debilidades</h2></div>"
             "<p>Solo las variables con dato. Las 34 restantes están pendientes y no castigan el índice.</p></div><div class='sw'>"
             f"<div class='card swc ok'><h3><i></i>Fortalezas</h3><div class='s'>Puntuación ≥ 80 %</div>{lst(lambda s: s >= .8, 'Ninguna variable alcanza 80 %.')}</div>"
             f"<div class='card swc mid'><h3><i></i>Se acerca</h3><div class='s'>Puntuación 50–80 %</div>{lst(lambda s: .5 <= s < .8, 'Ninguna.')}</div>"
             f"<div class='card swc bad'><h3><i></i>Debilidades</h3><div class='s'>Puntuación &lt; 50 %</div>{lst(lambda s: s < .5, 'Ninguna variable con dato cae bajo 50 %.')}</div>"
             "</div></section>")

    # VARIABLES
    H.append("<section class='blk' id='variables'><div class='sh'><div><div class='n'>05 · Propiedades</div><h2>Las 37 variables que puntúan</h2></div>"
             "<p>Estructura del Diccionario DAE. Pulsa SRC para ver método y fuente. La columna Peso es R × F.</p></div><div class='tools'>")
    for est in ("dentro", "cercano", "fuera", "sin dato"):
        H.append(f"<button type='button' class='tbtn' data-f='{est}'><span class='dot'></span><b>{pend_cnt[est]}</b>{fc.ESTADO_TXT[est]}</button>")
    H.append("<span class='sp'></span><input id='q' class='search' type='search' placeholder='Buscar variable…'>"
             "<button type='button' class='btn' id='openall'>Expandir todo</button><button type='button' class='btn' id='closeall'>Contraer</button>"
             "<button type='button' class='btn' id='clearf'>Limpiar</button></div>")
    grupos = [("Puntuación proxy", "Proxy · con dato QGIS"), ("Puntuación directa", "Directas"), ("Puntuación categórica", "Categóricas")]
    for k, (cat, titulo) in enumerate(grupos, 1):
        rs = [r for r in mat if r["Categoria"] == cat]
        con = sum(1 for r in rs if r["Variable"] in vals)
        sc = [vals[r["Variable"]][1] for r in rs if r["Variable"] in vals]
        media = sum(sc) / len(sc) if sc else None
        H.append(f"<div class='card acc{' open' if cat == 'Puntuación proxy' else ''}' id='sec-{k}'><button type='button' class='acch' aria-expanded='{str(cat == 'Puntuación proxy').lower()}'>"
                 f"<span class='ix'>{k:02d}</span><span class='t'>{e(titulo)}<small>{len(rs)} variables · {con} con dato</small></span>"
                 f"<span class='minis'>{'<span class=\"mini ok\">' + str(con) + ' con dato</span>' if con else ''}<span class='mini na'>{len(rs) - con} —</span></span>"
                 f"<span class='scw'><b>{(f'{media * 100:.0f}%' if media is not None else '—')}</b>{fc._bar(media)}</span><span class='chev'></span></button>"
                 "<div class='accb'><div><div class='vh'><span>Variable</span><span>Valor del sitio</span><span>Ventana Marte</span><span>Estado</span>"
                 "<span style='text-align:right'>Punt.</span><span style='text-align:center'>Peso</span><span></span></div>")
        for r in rs:
            if r["Variable"] in vals:
                v, s, nota = vals[r["Variable"]]
                H.append(fila_var(r, v, s, nota, POR_QUE.get(r["Variable"], "")))
            else:
                H.append(fila_var(r, "Pendiente", None, "", ""))
        H.append("</div></div></div>")
    H.append("</section>")

    # LÍMITES + CONTEXTO
    lim = [r for r in mat if r["Categoria"] == "Limitante"]
    H.append("<section class='blk rv' id='limites'><div class='sh'><div><div class='n'>06 · Límites</div><h2>Lo que la Tierra no puede imitar</h2></div>"
             "<p>No puntúan. Se declaran para que el límite del análogo quede explícito.</p></div><div class='grid3'>")
    for r in lim:
        sitio = ""
        if r["Variable"].startswith("Latitud"):
            sitio = f"Sitio: {f['lat']:.2f}° N (dentro del rango numérico, pero la latitud terrestre no reproduce la marciana)"
        H.append(f"<div class='card lim'><div class='t'><span>{e(r['Variable'])}</span></div><div class='mv'>{e(r['Referencia_Marte'])}</div>"
                 f"<div class='ml'>Valor en Marte</div>" + (f"<div class='sv'>{e(sitio)}</div>" if sitio else "") + f"<p>{e(r['Metodo_o_fuente'])}</p></div>")
    H.append("</div>")
    ctx = [r for r in mat if r["Categoria"] == "Contexto/misión"]
    H.append("<div class='sub2'><h3>Contexto, viabilidad y misión</h3><span>no puntúan</span></div><div class='card list'>"
             "<div class='lr hd'><span>Variable</span><span>Valor / criterio</span><span>Fuente</span></div>")
    for r in ctx:
        val = r["Referencia_Marte"]
        if r["Variable"].startswith("Extensión") and g:
            ok = g["extension_km"] >= 10
            val = f"{g['extension_km']:,.1f} km — {'cumple' if ok else 'no cumple'} el criterio ≥ 10 km"
        fuente = CAPAS_MX.get(r["Variable"], r["Metodo_o_fuente"])
        H.append(f"<div class='lr'><span class='k'>{e(r['Variable'])}</span><span class='m'>{e(val)}</span><span class='s'>{e(fuente)}</span></div>")
    H.append("</div></section>")

    # SPACEPORT
    la = abs(f["lat"]); s_lat = 1.0 if la <= 28 else max(0.0, 1 - (la - 28) / 10)
    crit = [("Latitud (valor absoluto)", f"{la:.2f}°", "0° a 28°", estado(s_lat), "Calculada del centroide del polígono"),
            ("Pendiente del terreno", "—", "≤ 5°", "sin dato", "CEM INEGI / SRTM"),
            ("Capacidad portante de la roca", "—", "> 300 kPa", "sin dato", "Ensayo de campo"),
            ("Resistencia de la roca", "—", "100–300 MPa", "sin dato", "SGM / ensayo"),
            ("Humedad relativa", "—", "50–70 %", "sin dato", "NASA POWER RH2M"),
            ("Días de lluvia al año", "—", "< 50–60", "sin dato", "Normales SMN"),
            ("Velocidad del viento", "—", "< 7–10 m/s", "sin dato", "NASA POWER WS10M"),
            ("Cizalladura del viento", "—", "< 0.05 kn/ft", "sin dato", "Radiosondeo / reanálisis"),
            ("Profundidad del nivel freático", "—", "≥ 3–5 m", "sin dato", "CONAGUA"),
            ("Presión hidrostática máxima", "—", "< 50–100 kPa", "sin dato", "CONAGUA / geotecnia"),
            ("Permeabilidad hidráulica", "—", "≥ 0.001 m/s", "sin dato", "Hidrogeología CONABIO"),
            ("Temperatura atmosférica", "—", "Moderada a cálida", "cualitativo", "Bioclimas CONABIO")]
    H.append("<section class='blk rv' id='spaceport'><div class='sh'><div><div class='n'>07 · Operaciones en Tierra</div><h2>Aptitud como spaceport terrestre</h2></div>"
             "<p>Índice aparte del DAE. Con 1 de 12 criterios evaluados no se calcula todavía.</p></div>"
             f"<div class='card sp-head'><div class='num'>1/12<small>criterios evaluados</small></div><div><p>Latitud {la:.2f}°: "
             + ("dentro del rango 0–28° que aprovecha la rotación terrestre." if s_lat >= .999 else f"{la - 28:.2f}° al norte del límite de 28° (puntuación {s_lat * 100:.0f} %).")
             + f"</p>{fc._bar(s_lat)}</div></div><div class='card list'>"
             "<div class='spr lr hd' style='display:grid'><span>Criterio</span><span>Valor</span><span>Umbral</span><span>Estado</span></div>")
    for a, b, c, est, src in crit:
        H.append(f"<div class='spr'><div class='k'><span>{e(a)}</span>{fc._tip(src, 'Fuente')}</div><div class='m'>{e(b)}</div><div class='c'>{e(c)}</div><div>{fc._pill(est)}</div></div>")
    H.append("</div></section>")

    # MÉTODO
    notas = [
        "Dominio marciano: zonas candidatas a exploración humana (latitud ±50°, elevación ≤ +2 km MOLA).",
        "Puntuación: 1 dentro de la ventana marciana; fuera decae linealmente: s = max(0, 1 − d/T).",
        "s_P = max(0, 1 − P/2); s_N = 1 si NDVI ≤ 0.1, si no max(0, 1 − (NDVI − 0.1)/0.3); s_T = max(0, 1 − (40 − T)/60).",
        "DAE* = 100 (0.30 s_P + 0.30 s_N + 0.40 s_T). El 40 % se asigna provisionalmente a la amplitud; no representa por sí sola el bloque restante.",
        "Cobertura = 60 % + 40 % × (6/541): fracción del peso metodológico con dato. No es parecido con Marte.",
        "Las variables sin medición no reciben cero: quedan como pendientes y solo bajan la cobertura.",
        "El índice QGIS (ind_mean) y el DAE se construyen con criterios distintos y no deben presentarse como equivalentes.",
        "Escala DAE-MX: nota = 1 + 9 × (P_QGIS + P_DAE)/200, con P = % de las 921 zonas que la zona supera. Alto ≥ 9; medio-alto 7–9; moderado 5–7; bajo < 5. Es una escala relativa al país.",
        "Conformidad: criterios de la hoja «Rangos de Variables» del equipo con dato disponible; cada uno suma su escala (5 o 4) si el valor cae en el rango.",
        "Lugares de referencia: documento de zonas del equipo; los 20 puntos caen dentro de su polígono. En la décima zona se corrigió −144.199 por −114.199.",
        f"Nota del cálculo: {f['nota']}",
        "Rásters de QGIS: EPSG:4326, ~0.0083° (~1 km). Estadística zonal recalculada por píxel; las medias coinciden con QGIS.",
    ]
    H.append("<section class='blk rv' id='metodo'><div class='sh'><div><div class='n'>Método</div><h2>Notas metodológicas</h2></div></div>"
             "<div class='card notes'><ul>" + "".join(f"<li>{e(x)}</li>" for x in notas) + "</ul></div></section></main>")
    H.append("<footer><div class='wrap'><div class='b'><i></i>Dynamic Analog Engine · Análogos México</div>"
             f"<p>FID {f['fid']} · Datos QGIS y Diccionario DAE · NASA Space Apps Challenge 2026</p></div></footer>"
             "<button id='up' type='button' aria-label='Volver arriba'></button><div id='tt' role='tooltip'></div>"
             f"<script>{fc.FICHA_JS}(function(){{document.querySelectorAll('.stack i[data-w]').forEach(function(x){{setTimeout(function(){{x.style.width=x.dataset.w+'%'}},300)}})}})();</script></body></html>")
    return "".join(H)


def scatter(fichas, w=900, h=440):
    x0, x1, y0, y1 = 5.5, 10.0, 62.0, 82.0
    L, Rr, T, B = 70, 30, 24, 54

    def X(v): return L + (v - x0) / (x1 - x0) * (w - L - Rr)

    def Y(v): return T + (y1 - v) / (y1 - y0) * (h - T - B)
    o = [f"<svg viewBox='0 0 {w} {h}' class='scat' role='img' aria-label='Nota nacional frente a DAE'>",
         f"<rect class='band' x='{X(9):.1f}' y='{T}' width='{X(10) - X(9):.1f}' height='{h - T - B}'/>"]
    for v in (6, 7, 8, 9, 10):
        o.append(f"<line class='gd' x1='{X(v):.1f}' x2='{X(v):.1f}' y1='{T}' y2='{h - B}'/><text x='{X(v):.1f}' y='{h - B + 18}' text-anchor='middle'>{v}</text>")
    for v in (65, 70, 75, 80):
        o.append(f"<line class='gd' x1='{L}' x2='{w - Rr}' y1='{Y(v):.1f}' y2='{Y(v):.1f}'/><text x='{L - 10}' y='{Y(v) + 4:.1f}' text-anchor='end'>{v} %</text>")
    o.append(f"<line class='ax' x1='{L}' x2='{w - Rr}' y1='{h - B}' y2='{h - B}'/><line class='ax' x1='{L}' x2='{L}' y1='{T}' y2='{h - B}'/>")
    o.append(f"<text x='{(L + w - Rr) / 2:.0f}' y='{h - 10}' text-anchor='middle'>Nota nacional DAE-MX (1–10) · franja verde = análogo alto</text>")
    o.append(f"<text x='16' y='{(T + h - B) / 2:.0f}' transform='rotate(-90 16 {(T + h - B) / 2:.0f})' text-anchor='middle'>DAE frente a Marte</text>")
    pts = sorted(fichas, key=lambda f: Y(f["dae"]))
    prev = -1e9
    for f in pts:
        cx, cy = X(min(f["nota"], 9.999)), Y(f["dae"])
        o.append(f"<circle class='pt-c' cx='{cx:.1f}' cy='{cy:.1f}' r='6'><title>{e(f['nombre'])}: nota {f['nota']:.2f}, DAE {f['dae']:.2f} %</title></circle>")
        if f["nota"] >= 9:
            ty = max(cy + 4, prev + 17); prev = ty
            tx = X(9) - 24
            o.append(f"<line x1='{tx + 4:.1f}' y1='{ty - 4:.1f}' x2='{cx - 7:.1f}' y2='{cy:.1f}' stroke='rgba(255,217,196,.35)' stroke-width='1'/>")
            o.append(f"<text class='lbt' x='{tx:.1f}' y='{ty:.1f}' text-anchor='end'>{f['rank_nota']}. {e(f['nombre'])} · {f['nota']:.2f}</text>")
        else:
            o.append(f"<text class='lbt' x='{cx + 12:.1f}' y='{cy + 4:.1f}'>{f['rank_nota']}. {e(f['nombre'])} · {f['nota']:.2f}</text>")
    o.append("</svg>")
    return "".join(o)


def index_html(fichas, otras, mapa, N):
    por_dae = sorted(fichas, key=lambda f: f["rank_nota"])
    etiquetas = {f["fid"]: str(f["rank_nota"]) for f in fichas}
    jac = [f for f in fichas if f["rank_dae_nac"] > 100]
    H = ["<!DOCTYPE html><html lang='es'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>",
         "<title>DAE · Fichas científicas de análogos marcianos en México</title>",
         f"<style>{fc.FICHA_CSS}{EXTRA_CSS}</style></head><body>",
         "<header class='top'><div class='wrap'><div class='brand'><i></i>Dynamic Analog Engine</div><nav class='nav'>"
         "<a href='#ranking'>Ranking</a><a href='#mapa'>Mapa</a><a href='#hallazgos'>Hallazgos</a><a href='#metodo'>Método</a></nav></div></header>"]
    best = por_dae[0]
    n_alto = sum(1 for f in fichas if f["nota"] >= 9)
    H.append("<section class='hero'><div class='wrap'><div class='hgrid'><div><div class='kick'>NASA Space Apps 2026 · Análogos México</div>"
             "<h1>Fichas científicas de 10 análogos marcianos</h1>"
             f"<p class='desc'>Las 10 zonas mejor calificadas por QGIS entre {N} candidatas, evaluadas con el Diccionario DAE: analogía, cobertura y confianza por separado.</p>"
             f"<div class='chips'><span class='chip'>{N} zonas candidatas</span><span class='chip'>3 de 37 variables con dato</span><span class='chip'>Baja California · Sonora · BCS</span></div></div>"
             f"<div class='gauge'><svg viewBox='0 0 120 120'><defs><linearGradient id='gg' x1='0' y1='0' x2='1' y2='1'><stop offset='0' stop-color='#e0782f'/><stop offset='1' stop-color='#c1440e'/></linearGradient></defs>"
             f"<circle class='trk' cx='60' cy='60' r='54' fill='none' stroke-width='7'/><circle class='val' data-v='{n_alto / 10:.4f}' cx='60' cy='60' r='54' fill='none' stroke-width='7'/></svg>"
             f"<div class='ctr'><div class='big' data-count='{n_alto}'>{n_alto}</div><div class='lab'>de 10 · análogo alto</div><div class='sub'>escala nacional DAE-MX</div></div></div></div>"
             "<div class='kpis'>"
             f"<div class='kpi'><div class='l'>Nota más alta (1–10)</div><div class='v'>{best['nota']:.2f}</div><div class='d'>{e(best['nombre'])}</div></div>"
             f"<div class='kpi'><div class='l'>Mejor DAE frente a Marte</div><div class='v'>{max(f['dae'] for f in fichas):.2f}<small>%</small></div><div class='d'>{e(max(fichas, key=lambda f: f['dae'])['nombre'])}</div></div>"
             f"<div class='kpi'><div class='l'>Cobertura de datos</div><div class='v'>{fichas[0]['cobertura']:.2f}<small>%</small></div><div class='d'>Igual en las 10 zonas</div></div>"
             f"<div class='kpi'><div class='l'>Confianza</div><div class='v'>{fichas[0]['confianza']:.0f}<small>%</small></div><div class='d'>3 variables proxy</div></div>"
             "</div></div></section><main class='wrap'>")
    H.append("<section class='blk rv' id='ranking'><div class='sh'><div><div class='n'>01 · Escala 1–10</div><h2>Top 10 en la escala nacional DAE-MX</h2></div>"
             "<p>Nota 1–10 relativa a las 921 zonas candidatas de México. Al lado, el DAE frente a Marte y la conformidad con los rangos del equipo.</p></div><div class='card list'>"
             "<div class='rk hd'><span>#</span><span>Zona</span><span></span><span>Nota</span><span>DAE</span><span>Rangos</span></div>")
    for f in por_dae:
        H.append(f"<a class='rk' href='fichas/{f['fid']}.html'><span class='n'>{f['rank_nota']:02d}</span><span class='t'><b>{e(f['nombre'])}</b><span>{e(f['clase_nac'])} · #{f['rank_qgis']} QGIS · #{f['rank_dae_nac']} DAE nacional · {len(f['lugares'])} {'lugar' if len(f['lugares']) == 1 else 'lugares'} de referencia</span></span>"
                 f"<span class='bar-c'>{fc._bar((f['nota'] - 1) / 9)}</span><span class='v'>{f['nota']:.2f}</span><span class='m'>{f['dae']:.1f} %</span><span class='m'>{f['conf']:.0f} %</span></a>")
    H.append("</div><div class='sub2'><h3>Dos escalas, una lectura</h3><span>nota nacional frente a DAE</span></div>"
             f"<div class='card mapcard'>{scatter(fichas)}</div>"
             "<p style='color:var(--mut);font-size:13px;margin-top:10px'>Derecha = más análogo dentro de México. Arriba = más cerca de las ventanas marcianas. "
             "Ninguna zona llega a 80 % de DAE mientras la amplitud térmica sea la única variable del bloque de 40 %.</p></section>")
    H.append("<section class='blk rv' id='mapa'><div class='sh'><div><div class='n'>02 · Mapa</div><h2>Dónde están</h2></div>"
             "<p>Números = posición en la escala nacional 1–10. Gris = el resto de las zonas candidatas.</p></div>"
             f"<div class='card mapcard'>{svg_mapa(mapa, etiquetas=etiquetas, w=900, h=640, puntos=mapa['puntos'])}"
             f"<div class='mapcap'><span><i style='background:#c1440e'></i>Top 10</span><span><i style='background:rgba(255,255,255,.18)'></i>{N} zonas candidatas</span></div></div></section>")
    H.append("<section class='blk rv' id='hallazgos'><div class='sh'><div><div class='n'>03 · Hallazgos</div><h2>Lo que cambia al aplicar el DAE</h2></div></div>")
    if jac:
        j = jac[0]
        H.append(f"<div class='callout'><b>{e(j['nombre'])}</b> es #{j['rank_qgis']} en QGIS pero #{j['rank_dae_nac']} de {N} con el DAE: "
                 f"su NDVI ({j['Nd']:.3f}) y su precipitación ({j['P']:.2f} mm/día) son los más altos del top 10. Conviene revisarla antes de presentarla como análogo.</div>")
    H.append("<div class='callout' style='border-color:var(--ok);background:rgba(62,207,142,.06)'>Las 10 zonas tienen amplitud térmica de 9 a 19 °C, muy lejos de los 40 °C marcianos. "
             "Ese componente limita el DAE de todas a menos de 80 % mientras sea la única variable del bloque de 40 %.</div>")
    if otras:
        H.append("<div class='sub2'><h3>Zonas fuera del top 10 QGIS con DAE alto</h3><span>candidatas a revisar</span></div><div class='card list' style='overflow-x:auto'><table class='zt'>"
                 "<tr><th>FID</th><th>DAE</th><th>DAE nacional</th><th>Ranking QGIS</th><th>NDVI</th><th>Precip. (mm/día)</th><th>Amplitud (Δ°C)</th></tr>"
                 + "".join(f"<tr><td>{o['fid']}</td><td>{o['dae']:.2f}%</td><td>#{o['rank_dae']}</td><td>#{o['rank_qgis']}</td><td>{o['ndvi']:.3f}</td><td>{o['precip']:.3f}</td><td>{o['rango']:.2f}</td></tr>" for o in otras)
                 + "</table></div>")
    H.append("</section>")
    H.append("<section class='blk rv' id='metodo'><div class='sh'><div><div class='n'>Método</div><h2>Tres resultados, no uno</h2></div></div><div class='idx'>"
             "<div class='card ic'><div class='l'>Escala nacional 1–10</div><p>Nota = 1 + 9 × (P_QGIS + P_DAE)/200, con P = % de las 921 zonas que la zona supera. Alto ≥ 9 (≈ 11 % superior del país).</p></div>"
             "<div class='card ic'><div class='l'>Analogía frente a Marte</div><p>DAE* = 100 (0.30 s_P + 0.30 s_N + 0.40 s_T): qué tan cerca están las variables con dato de sus ventanas marcianas.</p></div>"
             "<div class='card ic'><div class='l'>Conformidad con rangos del equipo</div><p>Criterios de la hoja «Rangos de Variables» con dato: NDVI, precipitación anual, amplitud, latitud y extensión. Cada uno suma su escala (5 o 4).</p></div>"
             "<div class='card ic'><div class='l'>Cobertura de datos</div><p>Cuánto del peso metodológico tiene dato: 60 % + 40 % × 6/541 = 60.44 %. Las variables pendientes no reciben cero.</p></div>"
             "<div class='card ic'><div class='l'>Confianza</div><p>Calidad del dato: campo o satélite directo = 1, proxy = 0.7, valor regional = 0.5. Hoy las 3 variables son proxies.</p></div>"
             "</div></section></main>"
             "<footer><div class='wrap'><div class='b'><i></i>Dynamic Analog Engine · Análogos México</div><p>Datos QGIS del equipo · Diccionario de variables Tierra–Marte · NASA Space Apps Challenge 2026</p></div></footer>"
             f"<button id='up' type='button' aria-label='Volver arriba'></button><div id='tt' role='tooltip'></div><script>{fc.FICHA_JS}</script></body></html>")
    return "".join(H)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--resumen", type=Path, required=True)
    ap.add_argument("--matriz", type=Path, required=True)
    ap.add_argument("--qgis", type=Path, required=True)
    ap.add_argument("--salida", type=Path, default=Path("salida_dae"))
    a = ap.parse_args()
    fichas, otras, mapa, N = preparar(a)
    (a.salida / "fichas").mkdir(parents=True, exist_ok=True)
    nombres = {f["fid"]: f["nombre"] for f in fichas}
    for f in fichas:
        (a.salida / "fichas" / f"{f['fid']}.html").write_text(ficha_html(f, mapa, nombres), encoding="utf-8")
    (a.salida / "index.html").write_text(index_html(fichas, otras, mapa, N), encoding="utf-8")
    limpio = [{k: v for k, v in f.items() if k not in ("polys", "matriz")} for f in fichas]
    (a.salida / "fichas_dae.json").write_text(json.dumps({"zonas": limpio, "fuera_top10_dae_alto": otras}, ensure_ascii=False, indent=1, default=float), encoding="utf-8")
    for f in sorted(fichas, key=lambda f: -f["dae"]):
        print(f"{f['fid']:>5} {f['nombre'][:38]:<38} DAE {f['dae']:6.2f} (csv {f['dae_csv']:6.2f})  nac #{f['rank_dae_nac']:<4} sens {f['sens_min']:.1f}-{f['sens_max']:.1f}  ext {f['geo']['extension_km']:.1f} km")


if __name__ == "__main__":
    main()
