#!/usr/bin/env python3
"""
Entregable: Fichas científicas de análogos marcianos en México.

Índices publicados
  IAM  Índice de Analogía Marciana (0–100 %): cercanía de las condiciones medidas a los rangos de referencia de Marte.
  ENA  Escala Nacional de Analogía (1–10): posición de la zona frente a las 921 zonas candidatas del país.
  IIG  Índice de Idoneidad Geoespacial (0–1): índice zonal del análisis SIG (QGIS).

Uso
  python fichas_iam.py --resumen resumen.csv --matriz matriz.csv --qgis carpeta/qgis --salida sitio
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI))
sys.path.insert(0, str(AQUI.parent / "fichas"))
import ficha_cientifica as fc  # noqa: E402
import fichas_dae as fd        # noqa: E402

e = fc.e

# --------------------------------------------------------------------------------------------
# Catálogo de variables con dato (explicación pública)
# --------------------------------------------------------------------------------------------
VARS = {
    "precip": dict(nombre="Precipitación media", unidad="mm/día", fuente="NASA GPM IMERG", res="0.1° (~10 km), remuestreado a ~1 km",
                   que="Lluvia promedio diaria sobre la zona.",
                   marte="0 mm/día. En Marte no llueve y el agua líquida no es estable en la superficie.",
                   porque="Un sitio hiperárido conserva suelos secos, sales y superficies parecidas a las marcianas, y permite probar equipos sin humedad.",
                   ventana="0 mm/día", puntua=True, peso="30 %"),
    "ndvi": dict(nombre="Índice de vegetación (NDVI)", unidad="adimensional", fuente="NASA MODIS MOD13A3", res="1 km",
                 que="Mide el verdor de las plantas. Valores de 0.1 o menos indican roca, arena o suelo desnudo.",
                 marte="0 a 0.1. Marte no tiene vegetación.",
                 porque="Las pruebas de rovers, aterrizaje y hábitats necesitan superficies desnudas como las de Marte.",
                 ventana="0 a 0.1", puntua=True, peso="30 %"),
    "rango": dict(nombre="Amplitud térmica diaria", unidad="°C", fuente="NASA MODIS MYD11A2 (temperatura de superficie, día − noche)", res="1 km",
                  que="Diferencia entre la temperatura de la superficie de día y de noche.",
                  marte="40 a 100 °C de oscilación diaria.",
                  porque="Los ciclos térmicos fuertes desgastan materiales, sellos, baterías y electrónica.",
                  ventana="40 a 100 °C", puntua=True, peso="40 %"),
    "polvo": dict(nombre="Polvo atmosférico", unidad="índice relativo", fuente="NASA MERRA-2", res="0.5° × 0.625° (~50 km), remuestreado",
                  que="Concentración de polvo en el aire cerca de la superficie.",
                  marte="Polvo permanente, con tormentas locales, regionales y globales.",
                  porque="Permite probar la limpieza de paneles solares y la protección de sellos y óptica.",
                  ventana="Contexto (no puntúa)", puntua=False, peso="no forma parte del cálculo; se reporta como contexto"),
}

# Áreas naturales protegidas (fuentes públicas; confirmar traslape con la capa federal de CONANP)
ANP = {
    785: "Colinda con la Reserva de la Biosfera Alto Golfo de California y Delta del Río Colorado (Mexicali y San Luis Río Colorado). Hay que confirmar el traslape y los permisos con CONANP.",
    105: "Isla Montague está en el área de la Reserva de la Biosfera Alto Golfo de California y Delta del Río Colorado y entre las islas del Golfo de California reconocidas por la UNESCO. El trabajo de campo requiere permisos de CONANP.",
    741: "Está en la franja de 28° a 30° N de Baja California que cubre el Área de Protección de Flora y Fauna Valle de los Cirios. Probablemente requiere permisos de CONANP.",
    781: "Está en la franja de 28° a 30° N de Baja California que cubre el Área de Protección de Flora y Fauna Valle de los Cirios. Probablemente requiere permisos de CONANP.",
    1011: "Está en la franja de 28° a 30° N de Baja California que cubre el Área de Protección de Flora y Fauna Valle de los Cirios. Probablemente requiere permisos de CONANP.",
    743: "Está cerca del límite noroeste del Área de Protección de Flora y Fauna Valle de los Cirios. Hay que confirmar el traslape con CONANP.",
    1397: "Abarca la franja del Valle de los Cirios (Baja California) y el municipio de Mulegé, donde está la Reserva de la Biosfera El Vizcaíno. Probable traslape con ambas áreas.",
    1197: "Bahía Tortugas es una comunidad dentro de la Reserva de la Biosfera El Vizcaíno. El trabajo de campo requiere permisos de CONANP.",
    1218: "Está en el municipio de Mulegé, que en su mayor parte pertenece a la Reserva de la Biosfera El Vizcaíno. Hay que confirmar el traslape con CONANP.",
}
ISLA = {105: "Es una isla sedimentaria en la desembocadura del río Colorado: el acceso es solo por mar y el terreno está influido por las mareas."}


def lcf(t):
    return t if t[:2].isupper() else t[:1].lower() + t[1:]


def clase_iam(i):
    return "alto" if i >= 80 else "moderado" if i >= 60 else "parcial" if i >= 40 else "débil"


# --------------------------------------------------------------------------------------------
# Cálculos adicionales
# --------------------------------------------------------------------------------------------
def enriquecer(fichas, args):
    nac = fd.leer_csv(args.qgis / "ranking_zonas.csv")
    iam_nac = np.array([fd.dae(float(z["lluvia_mean"]), float(z["ndvi_mean"]), float(z["rango_mean"])) for z in nac])
    for f in fichas:
        zs = f["zonal"]
        P3 = [zs["precip"]["p5"], f["P"], zs["precip"]["p95"]]
        N3 = [zs["ndvi"]["p5"], f["Nd"], zs["ndvi"]["p95"]]
        T3 = [zs["rango"]["p5"], f["T"], zs["rango"]["p95"]]
        W = (0.27, 0.30, 0.33)
        espacial = [fd.dae(p, n, t) for p, n, t in itertools.product(P3, N3, T3)]
        pesos = [fd.dae(f["P"], f["Nd"], f["T"], wp, wn) for wp, wn in itertools.product(W, W)]
        ambas = [fd.dae(p, n, t, wp, wn) for p, n, t, wp, wn in itertools.product(P3, N3, T3, W, W)]
        f["int"] = {"espacial": (min(espacial), max(espacial)), "pesos": (min(pesos), max(pesos)), "ambas": (min(ambas), max(ambas))}
        lo, hi = f["int"]["ambas"]

        def nota_de(iam):
            return 1 + 9 * (f["P_qgis"] + 100 * float(np.mean(iam_nac < iam))) / 200
        f["ena_int"] = (nota_de(lo), nota_de(hi))
        f["pct_top"] = (f["P_qgis"] + f["P_dae"]) / 2
        # fortalezas por variable (0–1) y percentil nacional "a favor de Marte"
        n = f["nac"]
        f["vars"] = [
            {"k": "ndvi", "valor": f"{f['Nd']:.3f}", "s": f["sN"], "pct": n["ndvi"]["pct_mas_verdes"],
             "pct_txt": f"Menos vegetación que el {n['ndvi']['pct_mas_verdes']:.0f} % de las zonas del país"},
            {"k": "precip", "valor": f"{f['P']:.2f} mm/día", "s": f["sP"], "pct": n["precip"]["pct_mas_humedas"],
             "pct_txt": f"Más seca que el {n['precip']['pct_mas_humedas']:.0f} % de las zonas del país"},
            {"k": "rango", "valor": f"{f['T']:.1f} °C", "s": f["sT"], "pct": n["rango"]["pct_menor"],
             "pct_txt": f"Mayor amplitud que el {n['rango']['pct_menor']:.0f} % de las zonas del país"},
            {"k": "polvo", "valor": f"{f['Pv'] * 1e8:.2f} × 10⁻⁸", "s": None, "pct": n["polvo"]["pct_menor"],
             "pct_txt": f"Más polvo que el {n['polvo']['pct_menor']:.0f} % de las zonas del país"},
        ]
        f["vars"].sort(key=lambda v: (-(v["s"] if v["s"] is not None else -1), -v["pct"]))
        f["favor"], f["limit"] = puntos(f)
    return iam_nac


def puntos(f):
    n = f["nac"]; fav, lim = [], []
    pct_top = 100 - f["pct_top"]
    if f["nota"] >= 9:
        fav.append(("Entre las zonas más análogas del país", f"En la Escala Nacional supera al {f['pct_top']:.1f} % de las 921 zonas candidatas de México."))
    if f["Nd"] <= 0.1:
        fav.append(("Superficie prácticamente desnuda", f"NDVI de {f['Nd']:.3f}, dentro del rango de referencia de Marte (0.1 o menos)."))
    elif f["Nd"] <= 0.13:
        fav.append(("Vegetación muy escasa", f"NDVI de {f['Nd']:.3f}, muy cerca del rango de referencia de Marte (0.1 o menos)."))
    if n["precip"]["pct_mas_humedas"] >= 90:
        fav.append(("Hiperaridez", f"{f['P']:.2f} mm/día de lluvia: más seca que el {n['precip']['pct_mas_humedas']:.0f} % de las zonas candidatas."))
    if n["polvo"]["pct_menor"] >= 80:
        fav.append(("Polvo atmosférico abundante", f"Más polvo que el {n['polvo']['pct_menor']:.0f} % de las zonas: buen escenario para probar la mitigación de polvo en paneles, sellos y óptica."))
    if n["rango"]["pct_menor"] >= 75:
        fav.append(("Oscilación térmica alta para México", f"{f['T']:.1f} °C diarios, mayor que en el {n['rango']['pct_menor']:.0f} % de las zonas."))
    if f["area_km2"] >= 1000:
        fav.append(("Superficie amplia", f"{f['area_km2']:,.0f} km² para recorridos largos de rovers y para desplegar infraestructura."))
    elif f["geo"]["extension_km"] >= 10:
        fav.append(("Extensión suficiente para aterrizajes", f"{f['geo']['extension_km']:.1f} km de extensión máxima, que alcanza para una elipse de aterrizaje de hasta 10 km."))
    if abs(f["lat"]) <= 28:
        fav.append(("Latitud favorable para lanzamientos", f"{f['lat']:.2f}° N, dentro del rango de 0 a 28° que aprovecha mejor la rotación terrestre."))
    if f["lugares"]:
        nom = ", ".join(L["nombre"].split(" (")[0] for L in f["lugares"][:3])
        fav.append(("Localidades de referencia dentro de la zona", f"{nom}: puntos de acceso para el trabajo de campo."))
    # limitaciones
    lim.append(("Oscilación térmica menor que la de Marte", f"{f['T']:.1f} °C frente a un mínimo de 40 °C en Marte. Es la mayor diferencia con Marte y resta {40 - f['aportes']['Amplitud térmica']:.1f} puntos al IAM."))
    if f["Nd"] > 0.13:
        lim.append(("Algo de vegetación", f"NDVI de {f['Nd']:.3f}, por encima del rango de referencia de Marte (0.1 o menos)."))
    if n["precip"]["pct_mas_humedas"] < 50:
        lim.append(("Más lluvia que el promedio de las candidatas", f"{f['P']:.2f} mm/día: más húmeda que la mitad de las zonas del país."))
    if f["area_km2"] < 100:
        lim.append(("Superficie reducida", f"{f['area_km2']:.0f} km², poco espacio para operaciones de gran escala o recorridos largos."))
    if f["fid"] in ISLA:
        lim.append(("Acceso marítimo", ISLA[f["fid"]]))
    if f["fid"] in ANP:
        lim.append(("Área natural protegida", ANP[f["fid"]]))
    if f["rank_dae_nac"] > 100:
        lim.append(("Analogía con Marte baja para el país", f"Su IAM la ubica en el lugar {f['rank_dae_nac']} de {f['N']}: la vegetación y la lluvia la alejan de Marte."))
    lim.append(("Variables por medir", "34 de las 37 variables comparables (geología, suelo, mineralogía, pendiente, cuevas) aún no tienen dato y requieren capas oficiales o campañas de campo."))
    return fav, lim


# --------------------------------------------------------------------------------------------
# Estilos y scripts adicionales
# --------------------------------------------------------------------------------------------
CSS2 = r"""
.brand .logo{height:32px;width:auto;margin-right:4px}footer .b .logo{height:22px;width:auto;margin-right:2px}
.imap .z-pl{fill:#fff;stroke:#c1440e;stroke-width:1.2;cursor:help}
.lead-x{font-size:18px;color:var(--ink2);max-width:760px;line-height:1.6;margin:0}
.explain{display:grid;grid-template-columns:repeat(3,1fr);gap:14px;margin-top:22px}
.step{padding:22px 24px;position:relative}.step .k{font:700 34px var(--disp);color:var(--burn);line-height:1}
.step h4{font:600 16px var(--disp);margin:10px 0 6px}.step p{font-size:13.5px;color:var(--ink2);margin:0;line-height:1.55}
.defs{display:grid;grid-template-columns:repeat(3,1fr);gap:14px}
.def{padding:24px}.def .ac{font:700 13px var(--mono);letter-spacing:.14em;color:var(--burn)}.def h4{font:700 20px var(--disp);margin:6px 0 10px}
.def p{font-size:13.5px;color:var(--ink2);margin:0 0 8px;line-height:1.55}.def .rg{font:500 12px var(--mono);color:var(--mut);border-top:1px solid var(--line);padding-top:10px;margin-top:12px}
.fv{display:grid;grid-template-columns:1fr 1fr;gap:16px}
.fl{padding:26px}.fl h3{display:flex;align-items:center;gap:10px;font:600 19px var(--disp);margin:0 0 16px}
.fl h3 i{width:12px;height:12px;border-radius:50%}.fl.ok h3 i{background:var(--ok);box-shadow:0 0 12px var(--ok)}.fl.bad h3 i{background:var(--bad);box-shadow:0 0 12px var(--bad)}
.fl .it{padding:14px 0;border-top:1px solid var(--line)}.fl .it:first-of-type{border-top:0}
.fl .it b{display:block;font:600 14.5px var(--disp);margin-bottom:3px}.fl .it span{font-size:13.5px;color:var(--ink2);line-height:1.55}
.vc-grid{display:grid;grid-template-columns:repeat(2,1fr);gap:14px}
.vc{padding:22px 24px;cursor:pointer;transition:.25s;position:relative}
.vc:hover{border-color:rgba(224,120,47,.5);transform:translateY(-2px)}
.vc .rk{position:absolute;top:18px;right:20px;font:700 12px var(--mono);color:var(--mut2);padding:0;border:0;display:block}
.vc .rk.top{color:var(--ok)}
.vc h4{font:600 17px var(--disp);margin:0 0 4px;padding-right:60px}.vc .src{display:block;font-size:12px;color:var(--mut);background:none;border:0;padding:0;margin:0}
.vc .nums{display:flex;align-items:baseline;gap:14px;margin:14px 0 6px}.vc .nums b{font:700 30px var(--disp)}.vc .nums span{font:500 12px var(--mono);color:var(--mut)}
.vc .two{display:grid;grid-template-columns:1fr 1fr;gap:12px;margin-top:12px}
.vc .two div{font:600 10.5px var(--mono);letter-spacing:.1em;text-transform:uppercase;color:var(--mut2)}
.vc .more{display:none;margin-top:16px;padding-top:14px;border-top:1px solid var(--line);font-size:13.5px;color:var(--ink2);line-height:1.6}
.vc .more p{margin:0 0 8px}.vc .more b{color:var(--ink)}
.vc.open .more{display:block}.vc .hint{font:500 11px var(--mono);color:var(--mut2);margin-top:12px}
.vc.open .hint{display:none}
.ci{padding:28px 30px}
.ci-tabs{display:flex;gap:8px;flex-wrap:wrap;margin-bottom:26px}
.ci-tabs button{font:500 13px var(--body);color:var(--ink2);background:var(--bg2);border:1px solid var(--line);padding:9px 14px;border-radius:10px;cursor:pointer;transition:.2s}
.ci-tabs button.on{color:var(--ink);border-color:var(--burn);background:rgba(193,68,14,.15)}
.ci-plot{position:relative;height:74px;margin:10px 0 4px}
.ci-plot .ax{position:absolute;left:0;right:0;top:36px;height:2px;background:var(--line2)}
.ci-plot .bandx{position:absolute;top:26px;height:22px;border-radius:7px;background:rgba(224,120,47,.28);border:1px solid rgba(224,120,47,.7);transition:left .5s,width .5s}
.ci-plot .ptx{position:absolute;top:22px;width:4px;height:30px;background:#fff;border-radius:2px;transform:translateX(-2px)}
.ci-plot .lab{position:absolute;top:0;font:600 12px var(--mono);transform:translateX(-50%);white-space:nowrap;transition:left .5s}
.ci-ax{display:flex;justify-content:space-between;font:500 11px var(--mono);color:var(--mut2)}
.ci-txt{font-size:14px;color:var(--ink2);margin:18px 0 0;line-height:1.6}
.ci-grid{display:grid;grid-template-columns:repeat(3,1fr);gap:14px;margin-top:14px}
.ci-grid .card{padding:20px 22px}.ci-grid .l{font:600 11px var(--mono);letter-spacing:.14em;text-transform:uppercase;color:var(--mut)}
.ci-grid .v{font:700 28px var(--disp);margin:8px 0 4px}.ci-grid p{font-size:13px;color:var(--ink2);margin:0}
.zmap{display:grid;grid-template-columns:1.4fr 1fr;gap:16px}
.zsvg{width:100%;height:auto;display:block}.zsvg .nb{fill:rgba(255,255,255,.07);stroke:rgba(255,255,255,.12);stroke-width:.4}
.zsvg .me{fill:rgba(193,68,14,.55);stroke:#ffd9c4;stroke-width:1}
.zsvg .lm circle{fill:#fff;stroke:#c1440e;stroke-width:2;cursor:pointer;transition:r .2s}.zsvg .lm:hover circle{r:8}
.zsvg .lm text{font:600 12px var(--body);fill:var(--ink);paint-order:stroke;stroke:#07080a;stroke-width:4px}
.zsvg .scl{stroke:var(--ink2);stroke-width:2}.zsvg .sct{font:500 11px var(--mono);fill:var(--ink2)}
.gloss{border-bottom:1px dotted var(--mut);cursor:help}
.srcs{display:grid;grid-template-columns:1.3fr 1fr;gap:16px}
.srcs h4{font:600 12px var(--mono);letter-spacing:.16em;text-transform:uppercase;color:var(--mut);margin:0 0 12px}
.srcs ul{list-style:none;margin:0;padding:0}.srcs li{padding:10px 0;border-top:1px solid var(--line);font-size:13.5px;color:var(--ink2);line-height:1.5}
.srcs li:first-child{border-top:0}.srcs li b{color:var(--ink);font-weight:600}.srcs a{color:var(--burn);text-decoration:none}.srcs a:hover{text-decoration:underline}
.sw-list .tag{display:inline-block;font:600 10.5px var(--mono);color:var(--burn);background:rgba(193,68,14,.12);padding:3px 8px;border-radius:6px;margin-right:8px}
.tbl{width:100%;border-collapse:collapse;font-size:13.5px}
.tbl th{font:600 10.5px var(--mono);letter-spacing:.12em;text-transform:uppercase;color:var(--mut2);text-align:left;padding:14px 16px;border-bottom:1px solid var(--line);background:var(--bg2);cursor:pointer;user-select:none;white-space:nowrap}
.tbl th:hover{color:var(--ink)}.tbl th.sorted{color:var(--burn)}.tbl th .ar{margin-left:6px;opacity:.6}
.tbl td{padding:14px 16px;border-bottom:1px solid var(--line);font-family:var(--mono);font-size:13px;white-space:nowrap}
.tbl td.nm{font-family:var(--body);font-weight:600;white-space:normal}.tbl td.nm a{text-decoration:none}.tbl td.nm a:hover{color:var(--burn)}
.tbl tr{transition:background .2s}.tbl tbody tr:hover{background:rgba(255,255,255,.03)}
.tbl .cl{font:600 11.5px var(--body);padding:4px 10px;border-radius:999px}
.cl.alto{color:var(--ok);background:var(--okbg)}.cl.medio{color:var(--burn);background:rgba(224,120,47,.14)}.cl.mod{color:var(--mid);background:var(--midbg)}.cl.bajo{color:var(--bad);background:var(--badbg)}
.sim{display:grid;grid-template-columns:320px 1fr;gap:20px}
.sim .ctl{padding:24px}.sim .ctl label{display:block;font:600 12px var(--mono);letter-spacing:.1em;text-transform:uppercase;color:var(--mut);margin:18px 0 8px}
.sim .ctl label:first-of-type{margin-top:0}.sim .ctl label b{float:right;color:var(--ink);font-size:14px}
.sim input[type=range]{width:100%;accent-color:#e0782f}
.sim .note{font-size:12.5px;color:var(--mut);margin-top:18px;line-height:1.55}
.sim .res{padding:20px 24px}.sbar{display:grid;grid-template-columns:28px 1fr 60px;gap:12px;align-items:center;padding:8px 0;transition:transform .4s}
.sbar .n{font:700 14px var(--disp);color:var(--burn)}.sbar .t{font-size:13.5px}.sbar .t small{display:block;color:var(--mut);font-size:11.5px}
.sbar .b{height:8px;background:rgba(255,255,255,.07);border-radius:9px;overflow:hidden;margin-top:6px}.sbar .b i{display:block;height:100%;background:linear-gradient(90deg,#c1440e,#e0782f);transition:width .4s}
.sbar .v{font:700 15px var(--disp);text-align:right}.sbar .d{font:500 11px var(--mono)}
.cmp{padding:26px}.cmp .sel{display:grid;grid-template-columns:1fr auto 1fr;gap:14px;align-items:center;margin-bottom:22px}
.cmp select{font:500 14px var(--body);color:var(--ink);background:var(--bg2);border:1px solid var(--line2);padding:10px 12px;border-radius:10px;width:100%}
.cmp .vs{font:700 14px var(--mono);color:var(--mut)}
.crow{display:grid;grid-template-columns:1fr 150px 1fr;gap:14px;align-items:center;padding:10px 0;border-top:1px solid var(--line)}
.crow .c{font:600 11px var(--mono);letter-spacing:.1em;text-transform:uppercase;color:var(--mut);text-align:center}
.crow .a,.crow .bb{display:flex;align-items:center;gap:10px}.crow .a{flex-direction:row-reverse}
.crow .bar2{flex:1;height:10px;background:rgba(255,255,255,.06);border-radius:9px;overflow:hidden;display:flex}
.crow .a .bar2{justify-content:flex-end}.crow .bar2 i{display:block;height:100%;transition:width .5s}
.crow .a .bar2 i{background:#e0782f}.crow .bb .bar2 i{background:#8fb3ff}
.crow .val{font:600 13px var(--mono);min-width:76px}.crow .a .val{text-align:right}
.crow .val.w{color:var(--ok)}
.tipbox{position:fixed;z-index:70;pointer-events:none;background:#1b1f27;border:1px solid var(--line2);border-radius:12px;padding:10px 14px;font-size:13px;color:var(--ink2);box-shadow:0 18px 40px -10px #000;opacity:0;transition:opacity .15s;max-width:280px}
.tipbox.on{opacity:1}.tipbox b{display:block;color:var(--ink);font:600 14px var(--disp);margin-bottom:3px}
.imap .z-top{cursor:pointer;transition:fill .2s}.imap .z-top:hover{fill:rgba(224,120,47,.85)}
.pending{columns:2;column-gap:28px}.pending .pg{break-inside:avoid;margin-bottom:18px}
.pending h5{font:600 12px var(--mono);letter-spacing:.14em;text-transform:uppercase;color:var(--burn);margin:0 0 8px}
.pending li{font-size:13.5px;color:var(--ink2);padding:6px 0;border-top:1px solid var(--line);list-style:none}.pending ul{margin:0;padding:0}
.pending li b{color:var(--ink);font-weight:600}.pending li span{display:block;font-size:12px;color:var(--mut)}
@media (max-width:980px){.explain,.defs,.fv,.vc-grid,.ci-grid,.zmap,.srcs,.sim{grid-template-columns:1fr}.pending{columns:1}
 .cmp .sel{grid-template-columns:1fr}.crow{grid-template-columns:1fr 78px 1fr;gap:8px}.crow .val{min-width:0;font-size:11.5px}.crow .bar2{min-width:36px}.tbl-wrap{overflow-x:auto}}
"""

JS_FICHA = r"""
(function(){
  const $=(s,c=document)=>c.querySelector(s),$$=(s,c=document)=>[...c.querySelectorAll(s)];
  $$('.vc').forEach(c=>c.addEventListener('click',()=>c.classList.toggle('open')));
  const ci=window.__CI;
  if(ci){
    const lo=ci.axis[0],hi=ci.axis[1],P=v=>((v-lo)/(hi-lo)*100).toFixed(2)+'%';
    const band=$('.bandx'),l1=$('#ci-lo'),l2=$('#ci-hi'),t=$('#ci-t');
    $('.ptx').style.left=P(ci.iam);
    function set(k){const [a,b]=ci[k];band.style.left=P(a);band.style.width=(((b-a)/(hi-lo))*100).toFixed(2)+'%';
      l1.style.left=P(a);l1.textContent=a.toFixed(1)+' %';l2.style.left=P(b);l2.textContent=b.toFixed(1)+' %';
      t.innerHTML=ci.txt[k];$$('.ci-tabs button').forEach(x=>x.classList.toggle('on',x.dataset.k===k));}
    $$('.ci-tabs button').forEach(b=>b.addEventListener('click',()=>set(b.dataset.k)));set('ambas');
  }
  const tip=$('.tipbox');
  $$('[data-tt]').forEach(el=>{el.addEventListener('mousemove',ev=>{tip.innerHTML=el.dataset.tt;tip.classList.add('on');
    tip.style.left=Math.min(ev.clientX+14,innerWidth-300)+'px';tip.style.top=(ev.clientY+14)+'px'});
    el.addEventListener('mouseleave',()=>tip.classList.remove('on'));});
})();
"""

JS_INDEX = r"""
(function(){
  const $=(s,c=document)=>c.querySelector(s),$$=(s,c=document)=>[...c.querySelectorAll(s)];
  const Z=window.__Z;
  // tabla ordenable
  const tb=$('#rank tbody');let dir={};
  $$('#rank th[data-k]').forEach(th=>th.addEventListener('click',()=>{const k=th.dataset.k;dir[k]=dir[k]===-1?1:-1;
    const rows=[...tb.rows].sort((a,b)=>{const x=a.dataset[k],y=b.dataset[k];const nx=parseFloat(x),ny=parseFloat(y);
      return (isNaN(nx)?x.localeCompare(y):nx-ny)*dir[k]});rows.forEach(r=>tb.appendChild(r));
    $$('#rank th').forEach(h=>{h.classList.toggle('sorted',h===th);const ar=h.querySelector('.ar');if(ar)ar.textContent=h===th?(dir[k]===1?'▲':'▼'):''});}));
  $$('.fchip').forEach(b=>b.addEventListener('click',()=>{const c=b.dataset.c,on=!b.classList.contains('on');
    $$('.fchip').forEach(x=>x.classList.remove('on'));if(on)b.classList.add('on');
    [...tb.rows].forEach(r=>r.style.display=(!on||r.dataset.cls===c)?'':'none');}));
  // simulador de pesos
  const sP=$('#wP'),sN=$('#wN');
  function sc(p,n,t,wp,wn){const wt=1-wp-wn;const a=Math.max(0,1-p/2),b=n<=0.1?1:Math.max(0,1-(n-0.1)/0.3),c=(t>=40&&t<=100)?1:Math.max(0,1-(t<40?40-t:t-100)/60);return 100*(wp*a+wn*b+wt*c)}
  function sim(){let wp=+sP.value/100,wn=+sN.value/100;if(wp+wn>0.95){wn=0.95-wp;sN.value=Math.round(wn*100)}
    $('#vP').textContent=Math.round(wp*100)+' %';$('#vN').textContent=Math.round(wn*100)+' %';$('#vT').textContent=Math.round((1-wp-wn)*100)+' %';
    const r=Z.map(z=>({z,v:sc(z.P,z.N,z.T,wp,wn)})).sort((a,b)=>b.v-a.v);const box=$('#simres');box.innerHTML='';
    r.forEach((o,i)=>{const d=(z=>z.rank_iam10)(o.z)-(i+1);const el=document.createElement('div');el.className='sbar';
      el.innerHTML=`<span class="n">${String(i+1).padStart(2,'0')}</span><span class="t">${o.z.nombre}<small>${d===0?'misma posición':(d>0?'▲ sube '+d:'▼ baja '+(-d))}</small><span class="b"><i style="width:${o.v.toFixed(1)}%"></i></span></span><span class="v">${o.v.toFixed(1)}%</span>`;box.appendChild(el)});}
  if(sP){sP.addEventListener('input',sim);sN.addEventListener('input',sim);$('#wreset').addEventListener('click',()=>{sP.value=30;sN.value=30;sim()});sim();}
  // comparador
  const A=$('#cA'),B=$('#cB');
  const MET=[['ena','Escala nacional',10,v=>v.toFixed(2),1],['iam','IAM',100,v=>v.toFixed(1)+' %',1],['iig','Idoneidad SIG',1,v=>v.toFixed(3),1],
    ['sP','Aridez (punt.)',1,v=>(v*100).toFixed(0)+' %',1],['sN','Sin vegetación (punt.)',1,v=>(v*100).toFixed(0)+' %',1],['sT','Amplitud (punt.)',1,v=>(v*100).toFixed(0)+' %',1],
    ['P','Lluvia mm/día',1.5,v=>v.toFixed(2),-1],['N','NDVI',0.25,v=>v.toFixed(3),-1],['T','Amplitud °C',25,v=>v.toFixed(1),1],['area','Área km²',null,v=>Math.round(v).toLocaleString('es-MX'),1]];
  function cmp(){const a=Z[A.value],b=Z[B.value];const box=$('#cmpres');box.innerHTML='';
    MET.forEach(([k,lab,mx,f,better])=>{const va=a[k],vb=b[k];const m=mx||Math.max(va,vb);const wa=better*(va-vb)>0,wb=better*(vb-va)>0;
      const el=document.createElement('div');el.className='crow';
      el.innerHTML=`<div class="a"><span class="val ${wa?'w':''}">${f(va)}</span><span class="bar2"><i style="width:${Math.min(100,va/m*100)}%"></i></span></div><div class="c">${lab}</div><div class="bb"><span class="bar2"><i style="width:${Math.min(100,vb/m*100)}%"></i></span><span class="val ${wb?'w':''}">${f(vb)}</span></div>`;box.appendChild(el)});}
  if(A){Z.forEach((z,i)=>{A.add(new Option(z.nombre,i));B.add(new Option(z.nombre,i))});A.value=0;B.value=Math.min(1,Z.length-1);A.onchange=cmp;B.onchange=cmp;cmp();}
  // mapa interactivo
  const tip=$('.tipbox');
  $$('[data-tt]').forEach(el=>{el.addEventListener('mousemove',ev=>{tip.innerHTML=el.dataset.tt;tip.classList.add('on');
    tip.style.left=Math.min(ev.clientX+14,innerWidth-300)+'px';tip.style.top=(ev.clientY+14)+'px'});
    el.addEventListener('mouseleave',()=>tip.classList.remove('on'));
    if(el.dataset.href)el.addEventListener('click',()=>location.href=el.dataset.href);});
})();
"""


# --------------------------------------------------------------------------------------------
# Fuentes y software (compartido)
# --------------------------------------------------------------------------------------------
FUENTES_DATOS = [
    ("NASA GPM IMERG", "Precipitación (0.1°).", "https://gpm.nasa.gov/data/imerg"),
    ("NASA MODIS MOD13A3", "Índice de vegetación NDVI mensual (1 km).", "https://appeears.earthdatacloud.nasa.gov"),
    ("NASA MODIS MYD11A2", "Temperatura de superficie día y noche (1 km).", "https://appeears.earthdatacloud.nasa.gov"),
    ("NASA MERRA-2", "Polvo atmosférico (reanálisis).", "https://giovanni.gsfc.nasa.gov"),
    ("NASA Earthdata", "Acceso a los productos satelitales.", "https://earthdata.nasa.gov"),
]
FUENTES_MARTE = [
    ("NASA Mars Fact Sheet", "Gravedad, presión, densidad atmosférica y duración del día marciano.", "https://nssdc.gsfc.nasa.gov"),
    ("NASA Mars Design Reference Architecture 5.0 (NASA-SP-2009-566)", "Capacidad portante, potencia, navegación y dosis de radiación.", "https://ntrs.nasa.gov"),
    ("NASA, Zonas de Exploración para misiones humanas (2016)", "Dominio de referencia: latitud ±50°, elevación ≤ +2 km.", "https://ntrs.nasa.gov/api/citations/20160002370/downloads/20160002370.pdf"),
    ("SWIM – Subsurface Water Ice Mapping (NASA/PSI)", "Hielo somero y fracción de hielo.", "https://swim.psi.edu/"),
    ("Mars Climate Database (LMD/ESA/CNES)", "Temperatura, viento, polvo e irradiancia en Marte.", "https://www-mars.lmd.jussieu.fr/mars/access.html"),
    ("REMS, Curiosity (PDS Atmospheres)", "Mediciones meteorológicas en superficie.", "https://pdsatmospheres.nmsu.edu/data_and_services/atmospheres_data/MARS/curiosity/rems.html"),
    ("Mars Orbital Data Explorer (PDS Geosciences)", "MOLA, THEMIS, TES, CRISM, HiRISE.", "https://ode.rsl.wustl.edu/mars"),
]
FUENTES_MX = [
    ("Geoportal CONABIO", "Edafología, climas, bioclimas, uso de suelo y regionalización.", "http://www.conabio.gob.mx/informacion/gis/"),
    ("CONANP – Áreas naturales protegidas federales", "Contexto de permisos (vía CONABIO).", "http://geoportal.conabio.gob.mx/metadatos/doc/html/anpenero2024gw.html"),
    ("INEGI – Geología", "Litología y estructuras.", "https://www.inegi.org.mx/temas/geologia/"),
    ("Servicio Geológico Mexicano – GeoInfoMex", "Geología 1:250 000.", "https://www.sgm.gob.mx/GeoInfoMexGobMx/"),
    ("Reserva de la Biosfera Alto Golfo de California (DOF)", "Decreto del área protegida.", "https://sidof.segob.gob.mx/notas/docFuente/4745874"),
    ("Valle de los Cirios (Protected Planet)", "Área de Protección de Flora y Fauna.", "https://protectedplanet.net/301939"),
]
SOFTWARE = [
    ("QGIS 3.40", "Delimitación de zonas, álgebra de mapas, estadística zonal e Índice de Idoneidad Geoespacial."),
    ("Python 3 · NumPy · Matplotlib · tifffile", "Lectura de rásters y GeoPackage, estadística por píxel, índices, intervalos y generación de las fichas."),
    ("NASA AppEEARS · Giovanni · Earthdata Search", "Descarga y preparación de productos satelitales."),
    ("Google Earth", "Inspección visual de las zonas."),
    ("HTML · CSS · JavaScript", "Fichas interactivas sin dependencias externas."),
]


def fuentes_html():
    def lst(items):
        return "<ul>" + "".join(f"<li><b><a href='{e(u)}' target='_blank' rel='noopener'>{e(n)}</a></b> — {e(d)}</li>" for n, d, u in items) + "</ul>"
    return ("<div class='srcs'><div class='card' style='padding:24px 26px'><h4>Datos de la Tierra (NASA)</h4>" + lst(FUENTES_DATOS) +
            "<h4 style='margin-top:22px'>Referencias de Marte</h4>" + lst(FUENTES_MARTE) +
            "<h4 style='margin-top:22px'>Fuentes oficiales mexicanas</h4>" + lst(FUENTES_MX) + "</div>"
            "<div class='card sw-list' style='padding:24px 26px'><h4>Software utilizado</h4><ul>" +
            "".join(f"<li><span class='tag'>SW</span><b>{e(n)}</b><br>{e(d)}</li>" for n, d in SOFTWARE) + "</ul></div></div>")


def head(titulo):
    return ("<!DOCTYPE html><html lang='es'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>"
            f"<title>Aztroteca · {e(titulo)}</title><style>{fc.FICHA_CSS}{fd.EXTRA_CSS}{CSS2}</style></head><body>")


def foot(js_extra=""):
    return (f"<footer><div class='wrap'><div class='b'><img src='{LOGO}' alt='Logo Aztroteca' class='logo'>Aztroteca · Análogos marcianos en México</div>"
            "<p>NASA Space Apps Challenge 2026 · Datos NASA y fuentes oficiales mexicanas</p></div></footer>"
            "<button id='up' type='button' aria-label='Volver arriba'></button><div id='tt' role='tooltip'></div><div class='tipbox'></div>"
            f"<script>{fc.FICHA_JS}{js_extra}</script></body></html>")


LOGO = "data:image/png;base64,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"
BRAND = f"<img src='{LOGO}' alt='Logo Aztroteca' class='logo'>Aztroteca"


def clase_css(c):
    return {"Análogo alto": "alto", "Análogo medio-alto": "medio", "Análogo moderado": "mod"}.get(c, "bajo")


# --------------------------------------------------------------------------------------------
# Mapa de la zona con lugares de referencia
# --------------------------------------------------------------------------------------------
def mapa_zona(f, mapa, w=640, h=460):
    g = f["polys"]
    xs = [x for p in g for r in p for x, _ in r] + [L["lon"] for L in f["lugares"]]
    ys = [y for p in g for r in p for _, y in r] + [L["lat"] for L in f["lugares"]]
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    mx, my = (x1 - x0) * .12 + .02, (y1 - y0) * .12 + .02
    x0, x1, y0, y1 = x0 - mx, x1 + mx, y0 - my, y1 + my
    k = math.cos(math.radians((y0 + y1) / 2))
    s = min(w / ((x1 - x0) * k), h / (y1 - y0))
    ox = (w - (x1 - x0) * k * s) / 2; oy = (h - (y1 - y0) * s) / 2

    def P(x, y):
        return ox + (x - x0) * k * s, oy + (y1 - y) * s

    def path(polys):
        return "".join("M" + "L".join("%.1f,%.1f" % P(x, y) for x, y in r) + "Z" for p in polys for r in p)
    nb = []
    for fid, gg in mapa["todas"].items():
        if fid == f["fid"]:
            continue
        bx = [x for p in gg for r in p for x, _ in r]; by = [y for p in gg for r in p for _, y in r]
        if not bx:
            continue
        if max(bx) >= x0 and min(bx) <= x1 and max(by) >= y0 and min(by) <= y1:
            nb.append(path(gg))
    o = [f"<svg viewBox='0 0 {w} {h}' class='zsvg' role='img' aria-label='Mapa de la zona'>",
         f"<path class='nb' d='{''.join(nb)}'/>",
         f"<path class='me' fill-rule='evenodd' d='{path(g)}' data-tt='<b>{e(f['nombre'])}</b>{f['area_km2']:,.0f} km²'/>"]
    for L in f["lugares"]:
        X, Y = P(L["lon"], L["lat"])
        nom = L["nombre"].split(" (")[0]
        anchor, dx = ("end", -11) if X > w * .65 else ("start", 11)
        o.append(f"<g class='lm' data-tt='<b>{e(nom)}</b>{L['lat']:.3f}°, {L['lon']:.3f}°'><circle cx='{X:.1f}' cy='{Y:.1f}' r='6'/>"
                 f"<text x='{X + dx:.1f}' y='{Y + 4:.1f}' text-anchor='{anchor}'>{e(nom)}</text></g>")
    # barra de escala
    km_px = s / 111.0
    paso = next(v for v in (1, 2, 5, 10, 20, 50, 100, 200) if v * km_px >= 60)
    o.append(f"<line class='scl' x1='20' y1='{h - 20}' x2='{20 + paso * km_px:.1f}' y2='{h - 20}'/><text class='sct' x='20' y='{h - 28}'>{paso} km</text>")
    o.append("</svg>")
    return "".join(o)


# --------------------------------------------------------------------------------------------
# Ficha
# --------------------------------------------------------------------------------------------
def ficha(f, mapa, N):
    i = f["dae"]; nota = f["nota"]; lo, hi = f["int"]["ambas"]; el, eh = f["ena_int"]
    clase = f["clase_nac"]
    H = [head(f"{f['nombre']} · Ficha científica de análogo marciano")]
    nav = [("resumen", "Resumen"), ("fortalezas", "Fortalezas"), ("balance", "A favor y limitaciones"), ("confiabilidad", "Confiabilidad"),
           ("ubicacion", "Ubicación"), ("datos", "Datos"), ("fuentes", "Fuentes")]
    H.append(f"<header class='top'><div class='wrap'><a class='brand' href='../index.html' style='text-decoration:none'>{BRAND}</a><nav class='nav'>"
             + "".join(f"<a href='#{a}'>{t}</a>" for a, t in nav) + "</nav></div></header>")
    # HERO
    H.append("<section class='hero' id='top'><div class='wrap'><div class='hgrid'><div>"
             f"<div class='kick'>Ficha científica · Zona {f['rank_nota']} de 10</div>"
             f"<h1>{e(f['nombre'])}</h1><p class='desc'>{e(f['descripcion'])}</p><div class='chips'>"
             f"<span class='chip'>{f['lat']:.3f}° N, {abs(f['lon']):.3f}° O</span><span class='chip'>{f['area_km2']:,.0f} km²</span>"
             f"<span class='chip'>Extensión {f['geo']['extension_km']:,.0f} km</span><span class='chip'>{len(f['lugares'])} {'localidad' if len(f['lugares']) == 1 else 'localidades'} de referencia</span></div>"
             f"<div class='cls'><i></i>{e(clase)} · Escala Nacional {nota:.2f}/10</div></div>"
             "<div class='gauge'><svg viewBox='0 0 120 120'><defs><linearGradient id='gg' x1='0' y1='0' x2='1' y2='1'>"
             "<stop offset='0' stop-color='#e0782f'/><stop offset='1' stop-color='#c1440e'/></linearGradient></defs>"
             "<circle class='trk' cx='60' cy='60' r='54' fill='none' stroke-width='7'/>"
             f"<circle class='val' data-v='{(nota - 1) / 9:.4f}' cx='60' cy='60' r='54' fill='none' stroke-width='7'/></svg>"
             f"<div class='ctr'><div class='big' data-count='{nota:.2f}' data-dec='2'>{nota:.2f}</div>"
             "<div class='lab'>Escala Nacional de Analogía</div><div class='sub'>de 1 a 10</div></div></div></div>")
    kp = [("Escala Nacional (ENA)", f"<span data-count='{nota:.2f}' data-dec='2'>{nota:.2f}</span><small>/10</small>", f"Intervalo {el:.2f}–{eh:.2f}", (nota - 1) / 9),
          ("Analogía con Marte (IAM)", f"<span data-count='{i:.1f}' data-dec='1'>{i:.1f}</span><small>%</small>", f"Intervalo {lo:.1f}–{hi:.1f} %", i / 100),
          ("Idoneidad SIG (IIG)", f"<span data-count='{f['qgis_index']:.3f}' data-dec='3'>{f['qgis_index']:.3f}</span>", f"Lugar {f['rank_qgis']} de {N}", (f["qgis_index"] - .73) / .08),
          ("Calidad de los datos", f"{f['confianza']:.0f}<small>%</small>", "Tres variables satelitales indirectas", f["confianza"] / 100),
          ("Cobertura de variables", f"{f['cobertura']:.0f}<small>%</small>", "Peso del índice con dato", f["cobertura"] / 100)]
    H.append("<div class='kpis k5'>" + "".join(
        f"<div class='kpi'><div class='l'>{t}</div><div class='v'>{v}</div><div class='d'>{e(d)}</div>{fc._bar(b, '')}</div>" for t, v, d, b in kp) + "</div></div></section><main class='wrap'>")

    # RESUMEN
    top = f["vars"][0]
    peor = min((v for v in f["vars"] if v["s"] is not None), key=lambda v: v["s"])
    resumen = [
        f"{f['nombre']} obtiene {nota:.2f} de 10 en la Escala Nacional de Analogía: supera al {f['pct_top']:.1f} % de las {N} zonas candidatas de México y se clasifica como {clase.lower()}.",
        f"Su Índice de Analogía Marciana es de {i:.1f} % (análogo {clase_iam(i)} frente a Marte), con un intervalo de confiabilidad de {lo:.1f} a {hi:.1f} %.",
        f"Su variable más fuerte es {lcf(VARS[top['k']]['nombre'])} ({top['valor']}). {top['pct_txt']}.",
        f"La mayor diferencia con Marte es {lcf(VARS[peor['k']]['nombre'])}: {peor['valor']} frente a {VARS[peor['k']]['ventana']} en Marte.",
    ]
    H.append("<section class='blk rv' id='resumen'><div class='sh'><div><div class='n'>01 · Resumen</div><h2>¿Qué tan parecida es a Marte?</h2></div></div>"
             "<div class='why'><div class='card story'>" + "".join(f"<p>{e(x)}</p>" for x in resumen) + "</div>"
             "<div class='card ic'><div class='l'>Escala Nacional de Analogía</div>"
             f"<div class='v'>{nota:.2f}<small>/ 10</small></div>" + fd.escala_html(nota) +
             f"<p>La escala compara la zona con las {N} zonas candidatas del país, combinando su analogía con Marte y su idoneidad geoespacial.</p></div></div></section>")

    # FORTALEZAS (variables)
    H.append("<section class='blk rv' id='fortalezas'><div class='sh'><div><div class='n'>02 · Fortalezas</div><h2>Variables más fuertes del sitio</h2></div>"
             "<p>Ordenadas de mayor a menor coincidencia con Marte. Toca una tarjeta para ver qué mide y por qué importa.</p></div><div class='vc-grid'>")
    for k, v in enumerate(f["vars"], 1):
        meta = VARS[v["k"]]
        s = v["s"]
        H.append(f"<div class='card vc'><span class='rk{' top' if k == 1 else ''}'>#{k}</span><h4>{e(meta['nombre'])}</h4><span class='src'>{e(meta['fuente'])}</span>"
                 f"<div class='nums'><b>{e(v['valor'])}</b><span>Marte: {e(meta['ventana'])}</span></div>"
                 f"<div class='two'><div>Coincidencia con Marte<br><b style='font:700 18px var(--disp);color:var(--ink);letter-spacing:0;text-transform:none'>{(f'{s * 100:.0f} %' if s is not None else 'Contexto')}</b>{fc._bar(s) if s is not None else ''}</div>"
                 f"<div>Posición en México<br><b style='font:700 18px var(--disp);color:var(--ink);letter-spacing:0;text-transform:none'>{v['pct']:.0f} %</b>{fc._bar(v['pct'] / 100, '')}</div></div>"
                 f"<div class='hint'>Ver explicación +</div><div class='more'><p><b>Qué mide:</b> {e(meta['que'])}</p><p><b>En Marte:</b> {e(meta['marte'])}</p>"
                 f"<p><b>Por qué importa:</b> {e(meta['porque'])}</p><p><b>En México:</b> {e(v['pct_txt'])}.</p><p><b>Resolución:</b> {e(meta['res'])} · <b>Peso en el IAM:</b> {e(meta['peso'])}</p></div></div>")
    H.append("</div></section>")

    # BALANCE
    H.append("<section class='blk rv' id='balance'><div class='sh'><div><div class='n'>03 · Evaluación del sitio</div><h2>Puntos a favor y limitaciones</h2></div></div><div class='fv'>"
             "<div class='card fl ok'><h3><i></i>Puntos a favor</h3>" + "".join(f"<div class='it'><b>{e(a)}</b><span>{e(b)}</span></div>" for a, b in f["favor"]) + "</div>"
             "<div class='card fl bad'><h3><i></i>Limitaciones</h3>" + "".join(f"<div class='it'><b>{e(a)}</b><span>{e(b)}</span></div>" for a, b in f["limit"]) + "</div></div></section>")

    # CONFIABILIDAD
    ax = (max(0, min(f["int"]["ambas"][0], i) - 6), min(100, max(f["int"]["ambas"][1], i) + 6))
    zs = f["zonal"]
    ci = {"iam": i, "axis": ax, "espacial": f["int"]["espacial"], "pesos": f["int"]["pesos"], "ambas": f["int"]["ambas"],
          "txt": {"espacial": f"Si se usan los valores bajos (P5) y altos (P95) de los {zs['rango']['n']} píxeles de la zona en lugar del promedio, el IAM varía entre {f['int']['espacial'][0]:.1f} y {f['int']['espacial'][1]:.1f} %.",
                  "pesos": f"Si el peso de la precipitación y el del NDVI cambian ±10 % (de 27 % a 33 %), el IAM varía entre {f['int']['pesos'][0]:.1f} y {f['int']['pesos'][1]:.1f} %.",
                  "ambas": f"Combinando la variabilidad dentro de la zona y los cambios de peso, el IAM se mantiene entre {lo:.1f} y {hi:.1f} %. La Escala Nacional correspondiente va de {el:.2f} a {eh:.2f}."}}
    H.append("<section class='blk rv' id='confiabilidad'><div class='sh'><div><div class='n'>04 · Confiabilidad</div><h2>Rangos de confiabilidad</h2></div>"
             "<p>Cuánto podría cambiar el resultado por la variación dentro de la zona y por la ponderación de las variables.</p></div>"
             "<div class='card ci'><div class='ci-tabs'><button type='button' data-k='ambas'>Combinado</button><button type='button' data-k='espacial'>Variación dentro de la zona</button>"
             "<button type='button' data-k='pesos'>Cambio de pesos ±10 %</button></div>"
             f"<div class='ci-plot'><span class='ax'></span><span class='bandx'></span><span class='ptx'></span><span class='lab' id='ci-lo'></span><span class='lab' id='ci-hi'></span></div>"
             f"<div class='ci-ax'><span>{ax[0]:.0f} %</span><span>IAM · valor central {i:.1f} %</span><span>{ax[1]:.0f} %</span></div><p class='ci-txt' id='ci-t'></p></div>"
             "<div class='ci-grid'>"
             f"<div class='card'><div class='l'>Calidad de los datos</div><div class='v'>{f['confianza']:.0f} %</div><p>Las tres variables vienen de satélite y son mediciones indirectas de lo que se mide en Marte. Una medición directa o de campo vale 100 %.</p></div>"
             f"<div class='card'><div class='l'>Cobertura de variables</div><div class='v'>{f['cobertura']:.0f} %</div><p>Parte del peso total del índice que ya tiene dato. Las variables sin medir no restan puntos.</p></div>"
             f"<div class='card'><div class='l'>Resolución</div><div class='v'>{zs['rango']['n']} px</div><p>Píxeles de ~1 km dentro de la zona. Con 4 o más píxeles el valor representa al sitio y no solo a la región.</p></div></div>"
             "<div class='sub2'><h3>Variación dentro de la zona</h3><span>estadística por píxel</span></div><div class='card list tbl-wrap'><table class='zt'>"
             "<tr><th>Variable</th><th>Píxeles</th><th>Promedio</th><th>Mediana</th><th>Bajo (P5)</th><th>Alto (P95)</th></tr>")
    for lab, kk, mult, nd in (("Precipitación (mm/día)", "precip", 1, 3), ("NDVI", "ndvi", 1, 3), ("Amplitud térmica (°C)", "rango", 1, 1), ("Polvo (× 10⁻⁸)", "polvo", 1e8, 2)):
        z = zs[kk]
        H.append(f"<tr><td>{lab}</td><td>{z['n']}</td>" + "".join(f"<td>{z[c] * mult:.{nd}f}</td>" for c in ("media", "mediana", "p5", "p95")) + "</tr>")
    H.append("</table></div></section>")

    # UBICACIÓN
    H.append("<section class='blk rv' id='ubicacion'><div class='sh'><div><div class='n'>05 · Ubicación</div><h2>Dónde está y cómo llegar</h2></div>"
             "<p>Pasa el cursor sobre los puntos para ver las coordenadas de cada localidad.</p></div>"
             f"<div class='zmap'><div class='card mapcard'>{mapa_zona(f, mapa)}</div><div><div class='card mapcard'>"
             + fd.svg_mapa(mapa, resaltar=f["fid"], puntos=mapa["puntos"]) +
             "<div class='mapcap'><span><i style='background:#c1440e'></i>Esta zona</span><span><i style='background:rgba(224,120,47,.5)'></i>Otras zonas del top 10</span></div></div>"
             "<div class='card list' style='margin-top:14px'>"
             + "".join(f"<div class='lr' style='grid-template-columns:1.4fr 1fr'><span class='k'>{e(L['nombre'].split(' (')[0])}</span><span class='m'>{L['lat']:.3f}°, {L['lon']:.3f}°</span></div>" for L in f["lugares"])
             + "</div></div></div></section>")

    # DATOS: Marte no reproducible + por medir
    lim = [r for r in f["matriz"] if r["Categoria"] == "Limitante"]
    pend = [r for r in f["matriz"] if r["Rol"] in ("Directa", "Categórica") or (r["Rol"] == "Proxy" and r["Variable"] not in ("Precipitación media", "NDVI", "Amplitud térmica diaria"))]
    grupos = {"Terreno y suelo": ["Pendiente", "Cobertura de rocas", "Tamaño", "Densidad", "Cohesión", "Ángulo", "Capacidad", "Hundimiento", "Textura"],
              "Química y agua": ["Percloratos", "Azufre", "Actividad", "pH", "Contenido de agua", "recurso hídrico", "Mineralogía", "Litología"],
              "Clima y energía": ["Temperatura", "Irradiancia", "viento", "Albedo", "potencia", "Opacidad", "Tormentas", "eólicos"],
              "Biología y subsuelo": ["Extremófilos", "Biofirmas", "cavidad", "tubo", "techo", "Estabilidad", "compresión"]}

    def grupo(v):
        for gname, keys in grupos.items():
            if any(kk.lower() in v.lower() for kk in keys):
                return gname
        return "Otras"
    por = {}
    for r in pend:
        por.setdefault(grupo(r["Variable"]), []).append(r)
    H.append("<section class='blk rv' id='datos'><div class='sh'><div><div class='n'>06 · Datos</div><h2>Lo que falta medir y lo que la Tierra no reproduce</h2></div>"
             "<p>Estas variables completarán el índice con capas oficiales o trabajo de campo. Las condiciones planetarias de Marte se declaran como límites de cualquier análogo terrestre.</p></div>"
             "<div class='card' style='padding:26px 28px'><div class='pending'>")
    for gname, rs in por.items():
        H.append(f"<div class='pg'><h5>{e(gname)}</h5><ul>" + "".join(
            f"<li><b>{e(r['Variable'])}</b><span>Marte: {e(r['Referencia_Marte'].split(';')[0].replace(' (provisional)', ''))} · {e(fd.CAPAS_MX.get(r['Variable'], r['Metodo_o_fuente']))}</span></li>" for r in rs) + "</ul></div>")
    H.append("</div></div><div class='sub2'><h3>Condiciones de Marte que ningún sitio terrestre reproduce</h3><span>límites de la analogía</span></div><div class='grid3'>")
    for r in lim:
        H.append(f"<div class='card lim'><div class='t'><span>{e(r['Variable'])}</span></div><div class='mv'>{e(r['Referencia_Marte'])}</div><div class='ml'>Valor en Marte</div></div>")
    H.append("</div></section>")

    H.append("<section class='blk rv' id='fuentes'><div class='sh'><div><div class='n'>07 · Fuentes</div><h2>Fuentes y software</h2></div></div>" + fuentes_html() + "</section></main>")
    H.append(foot("window.__CI=" + json.dumps(ci, ensure_ascii=False) + ";" + JS_FICHA))
    return "".join(H)


# --------------------------------------------------------------------------------------------
# Índice
# --------------------------------------------------------------------------------------------
def mapa_interactivo(mapa, fichas, w=900, h=640):
    svg = fd.svg_mapa(mapa, etiquetas={f["fid"]: str(f["rank_qgis"]) for f in fichas}, w=w, h=h, puntos=mapa["puntos"])
    # reemplaza el trazo conjunto del top por trazos individuales con tooltip y enlace
    todas = mapa["todas"]
    xs = [x for g in todas.values() for p in g for r in p for x, _ in r]; ys = [y for g in todas.values() for p in g for r in p for _, y in r]
    x0, x1, y0, y1 = min(xs) - .3, max(xs) + .3, min(ys) - .3, max(ys) + .3
    k = math.cos(math.radians((y0 + y1) / 2)); s = min(w / ((x1 - x0) * k), h / (y1 - y0))
    ox = (w - (x1 - x0) * k * s) / 2; oy = (h - (y1 - y0) * s) / 2
    piezas = []
    for f in fichas:
        d = "".join("M" + "L".join(f"{ox + (x - x0) * k * s:.1f},{oy + (y1 - y) * s:.1f}" for x, y in r) + "Z" for p in todas[f["fid"]] for r in p)
        tt = f"<b>Zona {f['rank_qgis']} · {e(f['nombre'])}</b>ENA {f['nota']:.2f}/10 (lugar {f['rank_nota']}) · IAM {f['dae']:.1f} %<br>Clic para abrir la ficha"
        piezas.append(f"<path class='z-top' d='{d}' data-tt=\"{e(tt)}\" data-href='fichas/{f['fid']}.html'/>")
    for f in fichas:
        for L in f["lugares"]:
            X = ox + (L["lon"] - x0) * k * s; Y = oy + (y1 - L["lat"]) * s
            tl = f"<b>{e(L['nombre'].split(' (')[0])}</b>Zona {f['rank_qgis']} · {L['lat']:.3f}°, {L['lon']:.3f}°"
            piezas.append(f"<circle class='z-pl' cx='{X:.1f}' cy='{Y:.1f}' r='3.2' data-tt=\"{e(tl)}\"/>")
    svg = svg.replace("<path class='z-top' d='", "<path class='z-top' style='display:none' d='", 1)
    import re as _re
    lbls = "".join(_re.findall(r"<g class='lbl'>.*?</g>", svg))
    svg = _re.sub(r"<g class='lbl'>.*?</g>", "", svg)
    return svg.replace("</svg>", "".join(piezas) + lbls + "</svg>", 1).replace("class='map'", "class='map imap'", 1)


def indice(fichas, mapa, N):
    fs = sorted(fichas, key=lambda f: f["rank_nota"])
    rank_iam = {f["fid"]: k for k, f in enumerate(sorted(fichas, key=lambda f: -f["dae"]), 1)}
    Z = [{"nombre": f["nombre"], "fid": f["fid"], "ena": f["nota"], "iam": f["dae"], "iig": f["qgis_index"], "sP": f["sP"], "sN": f["sN"], "sT": f["sT"],
          "P": f["P"], "N": f["Nd"], "T": f["T"], "area": f["area_km2"], "rank_iam10": rank_iam[f["fid"]]} for f in fs]
    n_alto = sum(1 for f in fichas if f["nota"] >= 9)
    H = [head("Análogos marcianos en México · Fichas científicas")]
    H.append(f"<header class='top'><div class='wrap'><div class='brand'>{BRAND}</div><nav class='nav'>"
             "<a href='#que-es'>Qué es</a><a href='#ranking'>Ranking</a><a href='#mapa'>Mapa</a><a href='#comparar'>Comparar</a><a href='#simulador'>Simulador</a>"
             "<a href='#datos'>Datos</a><a href='#limitaciones'>Limitaciones</a><a href='#fuentes'>Fuentes</a></nav></div></header>")
    H.append("<section class='hero'><div class='wrap'><div class='hgrid'><div><div class='kick'>Aztroteca · NASA Space Apps Challenge 2026</div>"
             "<h1>Análogos marcianos en México</h1>"
             f"<p class='desc'>Diez sitios del noroeste de México evaluados frente a las condiciones de Marte, a partir de datos satelitales de la NASA y fuentes oficiales mexicanas.</p>"
             f"<div class='chips'><span class='chip'>{N} zonas candidatas analizadas</span><span class='chip'>Baja California · Baja California Sur · Sonora</span></div></div>"
             "<div class='gauge'><svg viewBox='0 0 120 120'><defs><linearGradient id='gg' x1='0' y1='0' x2='1' y2='1'><stop offset='0' stop-color='#e0782f'/><stop offset='1' stop-color='#c1440e'/></linearGradient></defs>"
             f"<circle class='trk' cx='60' cy='60' r='54' fill='none' stroke-width='7'/><circle class='val' data-v='{n_alto / 10:.4f}' cx='60' cy='60' r='54' fill='none' stroke-width='7'/></svg>"
             f"<div class='ctr'><div class='big' data-count='{n_alto}'>{n_alto}</div><div class='lab'>de 10 · análogo alto</div><div class='sub'>Escala Nacional de Analogía</div></div></div></div>"
             "<div class='kpis'>"
             f"<div class='kpi'><div class='l'>Mejor sitio</div><div class='v' style='font-size:24px'>{e(fs[0]['nombre'])}</div><div class='d'>ENA {fs[0]['nota']:.2f}/10</div></div>"
             f"<div class='kpi'><div class='l'>IAM más alto</div><div class='v'>{max(f['dae'] for f in fichas):.1f}<small>%</small></div><div class='d'>{e(max(fichas, key=lambda f: f['dae'])['nombre'])}</div></div>"
             f"<div class='kpi'><div class='l'>Superficie evaluada</div><div class='v'>{sum(f['area_km2'] for f in fichas):,.0f}<small>km²</small></div><div class='d'>10 zonas</div></div>"
             f"<div class='kpi'><div class='l'>Localidades de referencia</div><div class='v'>{sum(len(f['lugares']) for f in fichas)}</div><div class='d'>todas dentro de su zona</div></div>"
             "</div></div></section><main class='wrap'>")
    # QUÉ ES
    H.append("<section class='blk rv' id='que-es'><div class='sh'><div><div class='n'>01 · Qué es</div><h2>Cómo medimos la analogía con Marte</h2></div></div>"
             "<p class='lead-x'>Un análogo terrestre es un lugar de la Tierra con condiciones parecidas a las de otro planeta. Sirve para probar rovers, hábitats, trajes e instrumentos antes de ir a Marte. "
             "Cada sitio se evalúa con tres indicadores complementarios.</p><div class='defs' style='margin-top:24px'>"
             "<div class='card def'><div class='ac'>IAM</div><h4>Índice de Analogía Marciana</h4><p>Mide, de 0 a 100 %, qué tan cerca están las condiciones medidas en el sitio de los rangos de referencia de Marte.</p>"
             "<p>Cada variable vale 100 % si cae dentro del rango marciano y baja de forma gradual conforme se aleja.</p><div class='rg'>Alto ≥ 80 % · Moderado 60–80 % · Parcial 40–60 %</div></div>"
             f"<div class='card def'><div class='ac'>ENA</div><h4>Escala Nacional de Analogía</h4><p>Ubica cada sitio, del 1 al 10, frente a las {N} zonas candidatas de México.</p>"
             "<p>Combina la posición del sitio en el IAM y en el Índice de Idoneidad Geoespacial: 10 significa que supera a todas las demás.</p><div class='rg'>Alto ≥ 9 · Medio-alto 7–9 · Moderado 5–7 · Bajo &lt; 5</div></div>"
             "<div class='card def'><div class='ac'>IIG</div><h4>Índice de Idoneidad Geoespacial</h4><p>Resultado del análisis espacial en QGIS que delimitó las zonas candidatas combinando lluvia, vegetación, temperatura y polvo.</p>"
             "<p>Valores más altos indican zonas más idóneas como análogo.</p><div class='rg'>Escala 0 a 1</div></div></div>"
             "<div class='explain'><div class='card step'><div class='k'>1</div><h4>Datos satelitales</h4><p>Precipitación, vegetación, temperatura de superficie y polvo de productos NASA a ~1 km.</p></div>"
             "<div class='card step'><div class='k'>2</div><h4>Comparación con Marte</h4><p>Cada variable se compara con su rango de referencia marciano para las zonas de exploración humana.</p></div>"
             "<div class='card step'><div class='k'>3</div><h4>Índices y confiabilidad</h4><p>Se calculan IAM, ENA e IIG, junto con un intervalo de confiabilidad para cada sitio.</p></div></div></section>")
    # RANKING
    H.append("<section class='blk rv' id='ranking'><div class='sh'><div><div class='n'>02 · Ranking</div><h2>Los 10 sitios</h2></div>"
             "<p>Toca un encabezado para ordenar. Filtra por clase o abre una ficha.</p></div>"
             "<div class='tools'><button type='button' class='tbtn fchip' data-c='alto' style='--c:var(--ok)'><span class='dot'></span>Análogo alto</button>"
             "<button type='button' class='tbtn fchip' data-c='mod' style='--c:var(--mid)'><span class='dot'></span>Análogo moderado</button></div>"
             "<div class='card list tbl-wrap'><table class='tbl' id='rank'><thead><tr>"
             "<th data-k='r'>#<span class='ar'></span></th><th data-k='n'>Sitio<span class='ar'></span></th><th data-k='ena'>ENA<span class='ar'></span></th>"
             "<th data-k='iam'>IAM<span class='ar'></span></th><th data-k='ci'>Intervalo IAM</th><th data-k='iig'>IIG<span class='ar'></span></th>"
             "<th data-k='ndvi'>NDVI<span class='ar'></span></th><th data-k='p'>Lluvia<span class='ar'></span></th><th data-k='a'>Área km²<span class='ar'></span></th><th>Clase</th></tr></thead><tbody>")
    for f in fs:
        lo, hi = f["int"]["ambas"]
        H.append(f"<tr data-r='{f['rank_nota']}' data-n='{e(f['nombre'])}' data-ena='{f['nota']:.4f}' data-iam='{f['dae']:.4f}' data-iig='{f['qgis_index']:.5f}' data-ndvi='{f['Nd']:.5f}' data-p='{f['P']:.5f}' data-a='{f['area_km2']:.2f}' data-cls='{clase_css(f['clase_nac'])}'>"
                 f"<td>{f['rank_nota']:02d}</td><td class='nm'><a href='fichas/{f['fid']}.html'>{e(f['nombre'])} →</a></td><td><b>{f['nota']:.2f}</b></td><td>{f['dae']:.1f} %</td><td>{lo:.1f}–{hi:.1f}</td>"
                 f"<td>{f['qgis_index']:.3f}</td><td>{f['Nd']:.3f}</td><td>{f['P']:.2f}</td><td>{f['area_km2']:,.0f}</td><td><span class='cl {clase_css(f['clase_nac'])}'>{e(f['clase_nac'].replace('Análogo ', ''))}</span></td></tr>")
    H.append("</tbody></table></div><div class='sub2'><h3>Dos lecturas de la analogía</h3><span>posición nacional frente a analogía con Marte</span></div>"
             f"<div class='card mapcard'>{fd.scatter(fichas).replace('DAE frente a Marte', 'IAM (%)').replace(', DAE ', ', IAM ').replace('frente a DAE', 'frente a IAM').replace(': nota ', ': ENA ').replace('Nota nacional DAE-MX (1–10)', 'Escala Nacional de Analogía (1–10)')}</div></section>")
    # MAPA
    H.append("<section class='blk rv' id='mapa'><div class='sh'><div><div class='n'>03 · Mapa</div><h2>Ubicación de los sitios</h2></div>"
             f"<p>Números = zonas 1 a 10 del documento de zonas (orden del Índice de Idoneidad Geoespacial). Puntos blancos = localidades de referencia. Pasa el cursor para ver detalles y haz clic para abrir la ficha.</p></div>"
             f"<div class='card mapcard'>{mapa_interactivo(mapa, fichas)}</div></section>")
    # COMPARADOR
    H.append("<section class='blk rv' id='comparar'><div class='sh'><div><div class='n'>04 · Comparar</div><h2>Compara dos sitios</h2></div>"
             "<p>En verde, el valor más favorable para un análogo marciano.</p></div>"
             "<div class='card cmp'><div class='sel'><select id='cA' aria-label='Sitio A'></select><span class='vs'>VS</span><select id='cB' aria-label='Sitio B'></select></div><div id='cmpres'></div></div></section>")
    # SIMULADOR
    H.append("<section class='blk rv' id='simulador'><div class='sh'><div><div class='n'>05 · Simulador</div><h2>¿Qué pasa si cambian los pesos?</h2></div>"
             "<p>Ajusta la importancia de cada variable y observa cómo cambia el IAM y el orden de los sitios.</p></div>"
             "<div class='sim'><div class='card ctl'><label>Precipitación <b id='vP'>30 %</b></label><input type='range' id='wP' min='0' max='80' value='30'>"
             "<label>Vegetación (NDVI) <b id='vN'>30 %</b></label><input type='range' id='wN' min='0' max='80' value='30'>"
             "<label>Amplitud térmica <b id='vT'>40 %</b></label><div class='bar'><i style='width:100%'></i></div>"
             "<p class='note'>La amplitud térmica recibe el peso restante. Valores de referencia: 30 % · 30 % · 40 %.</p>"
             "<button type='button' class='btn' id='wreset' style='margin-top:10px'>Restablecer</button></div><div class='card res' id='simres'></div></div></section>")
    # DATOS
    H.append("<section class='blk rv' id='datos'><div class='sh'><div><div class='n'>06 · Datos</div><h2>Variables evaluadas</h2></div>"
             "<p>Hoy hay datos de cuatro variables satelitales. Otras 34 variables comparables se integrarán con capas oficiales y trabajo de campo.</p></div><div class='vc-grid'>")
    for kk, meta in VARS.items():
        H.append(f"<div class='card vc open' style='cursor:default'><h4>{e(meta['nombre'])}</h4><span class='src'>{e(meta['fuente'])} · {e(meta['res'])}</span>"
                 f"<div class='more' style='display:block'><p><b>Qué mide:</b> {e(meta['que'])}</p><p><b>En Marte:</b> {e(meta['marte'])}</p><p><b>Por qué importa:</b> {e(meta['porque'])}</p>"
                 f"<p><b>Peso en el IAM:</b> {e(meta['peso'])}</p></div></div>")
    H.append("</div></section>")
    # LIMITACIONES GENERALES
    lims = ["Los resultados usan cuatro variables satelitales. Geología, suelo, mineralogía, pendiente y cuevas aún no tienen dato y pueden cambiar el resultado.",
            "La amplitud térmica de todos los sitios (9 a 19 °C) queda lejos de los 40 °C de Marte, así que ningún sitio alcanza 80 % de IAM con los datos actuales.",
            "Los valores satelitales son promedios de ~1 km por píxel; no sustituyen las mediciones de campo.",
            "La gravedad, la presión, la composición atmosférica, la radiación y la duración del día de Marte no se reproducen en ningún sitio terrestre.",
            "Varios sitios están en áreas naturales protegidas o colindan con ellas; cualquier actividad requiere confirmar permisos con CONANP."]
    H.append("<section class='blk rv' id='limitaciones'><div class='sh'><div><div class='n'>07 · Limitaciones</div><h2>Alcance de los resultados</h2></div></div>"
             "<div class='card notes'><ul>" + "".join(f"<li>{e(x)}</li>" for x in lims) + "</ul></div></section>")
    H.append("<section class='blk rv' id='fuentes'><div class='sh'><div><div class='n'>08 · Fuentes</div><h2>Fuentes y software</h2></div></div>" + fuentes_html() + "</section></main>")
    H.append(foot("window.__Z=" + json.dumps(Z, ensure_ascii=False) + ";" + JS_INDEX))
    return "".join(H)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--resumen", type=Path, required=True)
    ap.add_argument("--matriz", type=Path, required=True)
    ap.add_argument("--qgis", type=Path, required=True)
    ap.add_argument("--salida", type=Path, default=Path("sitio"))
    a = ap.parse_args()
    fichas, otras, mapa, N = fd.preparar(a)
    enriquecer(fichas, a)
    (a.salida / "fichas").mkdir(parents=True, exist_ok=True)
    for f in fichas:
        (a.salida / "fichas" / f"{f['fid']}.html").write_text(ficha(f, mapa, N), encoding="utf-8")
    (a.salida / "index.html").write_text(indice(fichas, mapa, N), encoding="utf-8")
    datos = [{k: v for k, v in f.items() if k not in ("polys", "matriz", "zonal")} for f in fichas]
    (a.salida / "datos.json").write_text(json.dumps(datos, ensure_ascii=False, indent=1, default=float), encoding="utf-8")
    for f in sorted(fichas, key=lambda f: f["rank_nota"]):
        print(f"{f['rank_nota']:>2} {f['nombre'][:36]:<36} ENA {f['nota']:.2f} [{f['ena_int'][0]:.2f}-{f['ena_int'][1]:.2f}]  IAM {f['dae']:.1f} [{f['int']['ambas'][0]:.1f}-{f['int']['ambas'][1]:.1f}]  favor {len(f['favor'])} lim {len(f['limit'])}")


if __name__ == "__main__":
    main()
