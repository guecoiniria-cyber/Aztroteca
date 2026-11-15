#!/usr/bin/env python3
"""
Generador de FICHAS CIENTÍFICAS de sitios análogos Tierra–Marte
=================================================================
Proyecto: Análogos en México / Dynamic Analog Engine (NASA Space Apps)

Entrada
  1. Capa de zonas ya limpia y poligonizada en QGIS (.geojson, .gpkg, .shp*, o .csv)
     Cada polígono = una zona candidata, con sus atributos (valores físicos por zona).
  2. Catálogo de variables (el CSV "Variables_Analogos_Tierra-Marte_corregido").
  3. mapeo_columnas.csv  → qué columna de la capa corresponde a cada variable
     (se genera automáticamente la primera vez con --crear-mapeo).
  4. criterios_spaceport.csv (opcional) → criterios terrestres de spaceport.

Salida (carpeta --salida)
  fichas/<id>.html   ficha científica por zona (lista para imprimir/PDF)
  fichas/<id>.md     misma ficha en Markdown
  fichas.json        todas las fichas en JSON (para el sitio en Vercel)
  ranking.csv        tabla resumen ordenada por índice de analogía
  zonas_resultado.geojson  polígonos + índices (para QGIS)
  zonas_resultado.kml      polígonos + ficha resumida en globo (Google Earth)
  index.html         índice navegable de todas las fichas

Uso
  python ficha_cientifica.py --crear-mapeo --catalogo Variables.csv
  python ficha_cientifica.py --zonas zonas.geojson --catalogo Variables.csv \
         --mapeo mapeo_columnas.csv --spaceport criterios_spaceport.csv --salida salida

Sin dependencias obligatorias (solo Python ≥3.9). Si geopandas está instalado se
usa para leer .shp / cualquier CRS; si no, .geojson y .gpkg se leen de forma nativa
(en EPSG:4326).
"""
from __future__ import annotations

import argparse
import csv
import html
import io
import json
import math
import re
import sqlite3
import struct
import sys
import unicodedata
from dataclasses import dataclass, field, asdict
from pathlib import Path

# ---------------------------------------------------------------------------
# Utilidades de texto
# ---------------------------------------------------------------------------

def norm(s: str) -> str:
    """minúsculas sin acentos, para comparar."""
    s = unicodedata.normalize("NFKD", str(s)).encode("ascii", "ignore").decode()
    return re.sub(r"\s+", " ", s.lower()).strip()


def slug(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", norm(s)).strip("_")[:60]


def leer_texto(path: Path) -> str:
    raw = path.read_bytes()
    for enc in ("utf-8-sig", "cp1252", "latin-1"):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            continue
    return raw.decode("utf-8", "replace")


def num(x):
    if x is None:
        return None
    if isinstance(x, (int, float)):
        return None if (isinstance(x, float) and math.isnan(x)) else float(x)
    s = str(x).strip().replace(",", ".")
    if s == "" or s.upper() in ("N/A", "NA", "NAN", "NULL", "NONE", "-"):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def fmt(v, nd=3):
    if v is None:
        return "—"
    if isinstance(v, str):
        return v
    if v == 0:
        return "0"
    a = abs(v)
    if a >= 1e5 or a < 1e-3:
        return f"{v:.3g}"
    if a >= 100:
        return f"{v:,.0f}".replace(",", " ")
    return f"{v:.{nd}g}" if a < 1 else f"{v:.{max(1, nd-1)}f}".rstrip("0").rstrip(".")


# Símbolos que Excel pierde al exportar CSV en Windows (se convierten en "?")
REPARAR = [
    ("relevancia ? 3", "relevancia ≥ 3"), ("relevancia ? 2", "relevancia ≤ 2"),
    ("rocas ?0.5", "rocas ≥0.5"), ("Opacidad atmosférica ?", "Opacidad atmosférica τ"),
    ("(SO?)", "(SO₃)"), ("[?°C | ?°F]", "[Δ°C | Δ°F]"), ("(? > 3)", "(τ > 3)"),
    ("día ? LST", "día − LST"), ("(?¹³C)", "(δ¹³C)"), ("(CO?)", "(CO₂)"),
    ("NDVI ? 0.1", "NDVI ≤ 0.1"), ("°C?°F", "°C→°F"), ("?°C??°F", "Δ°C→Δ°F"),
    ("aw < 0.6", "aw < 0.6"),
]


def reparar(s: str) -> str:
    for a, b in REPARAR:
        s = s.replace(a, b)
    s = re.sub(r"\s\?$", " ✓", s.strip())          # marca "fuente en la guía"
    return s


# ---------------------------------------------------------------------------
# Conocimiento por variable: columna sugerida, tolerancia, conversión e importancia
# (clave = palabra(s) que aparecen en el nombre de la variable del catálogo)
# ---------------------------------------------------------------------------
CONOCIMIENTO = [
    # clave, columna sugerida, tolerancia, nota de conversión, por qué importa
    ("litologia dominante", "litologia", None, "texto (clase de roca dominante)",
     "La roca basáltica es la corteza dominante de Marte; un sustrato equivalente permite probar perforación, abrasión, cimentación y calibrar espectrómetros (Raman, IR, DRX)."),
    ("mineralogia de alteracion", "mineralogia", None, "texto; EMIT L2B lista minerales",
     "Arcillas, sulfatos, sílice y hematita registran agua pasada en Marte: son blancos de astrobiología y de ISRU (agua de hidratación)."),
    ("tamano de grano", "grano_mm", None, "mm",
     "El tamaño de grano controla tracción, abrasión de sellos, filtración de polvo y la inercia térmica observada desde órbita."),
    ("densidad aparente", "densidad_kgm3", None, "kg/m³",
     "Interviene en capacidad portante, excavación y blindaje con regolito apilado."),
    ("cohesion", "cohesion_kpa", None, "kPa",
     "Junto con el ángulo de fricción define la resistencia al corte (Mohr-Coulomb) para taludes, zanjas y movilidad."),
    ("angulo de friccion", "friccion_deg", None, "°",
     "Parámetro clave de mecánica de suelos: estabilidad de bermas, excavación y tracción de ruedas."),
    ("capacidad portante", "portante_kpa", None, "kPa",
     "Determina si el terreno soporta módulos de hábitat, plataformas de aterrizaje y tanques sin cimentación profunda."),
    ("percloratos", "percloratos_wt", 1.0, "wt%",
     "Los percloratos son tóxicos para la tripulación, corrosivos y oxidantes de orgánicos; un análogo con sales oxidantes permite probar protocolos."),
    ("contenido de azufre", "azufre_wt", 5.0, "wt% SO₃",
     "El azufre es materia prima del concreto de azufre marciano (construcción ISRU) e indica ambientes evaporíticos."),
    ("elevacion del sitio", "elevacion_m", None, "m",
     "En Marte la elevación baja maximiza la columna de atmósfera para frenar el descenso (EDL)."),
    ("pendiente local", "pendiente_deg", 5.0, "° (SRTM 30 m, media de la zona)",
     "Pendientes ≤5° son requisito de aterrizaje y de construcción de hábitats y plataformas."),
    ("cobertura de rocas", "rocas_pct", 10.0, "% de área",
     "Las rocas ≥0.5 m pueden volcar un aterrizador o dañar ruedas; terreno plano a 30 m/píxel puede estar sembrado de bloques."),
    ("semieje mayor", "extension_km", None, "km (se calcula de la geometría si falta)",
     "La zona debe ser lo bastante grande para contener la elipse de aterrizaje y la huella del spaceport."),
    ("latitud del sitio", "latitud", None, "° (centroide)",
     "En Marte las latitudes medias combinan hielo somero accesible y energía solar suficiente."),
    ("pendiente transitable", "pendiente_max_deg", 10.0, "° (máximo o percentil 95)",
     "Pendientes de 15–20° comprometen tracción y estabilidad de rovers; define rutas de travesía."),
    ("hundimiento", "sinkage_cm", 10.0, "cm",
     "El hundimiento de rueda atrapó a Spirit en Marte; mide la transitabilidad real del suelo suelto."),
    ("precipitacion", "precip_mm_dia", 2.0, "mm/día (GPM IMERG mensual viene en mm/h → ×24)",
     "Marte no tiene lluvia líquida; la hiperaridez preserva el regolito seco, las sales y las biofirmas, y es el primer filtro del índice."),
    ("ndvi", "ndvi", 0.3, "adimensional (MOD13A3 viene ×10000 → factor 0.0001)",
     "Marte no tiene vegetación; NDVI ≤ 0.1 indica roca, arena o costra desnuda, condición necesaria para un análogo de superficie."),
    ("temperatura minima", "lst_noche_c", None, "°C (MOD11A2 viene en K×50 → factor 0.02, offset −273.15)",
     "El frío nocturno extremo gobierna el diseño térmico de baterías, electrónica y sellos; la Tierra solo se acerca."),
    ("temperatura maxima", "lst_dia_c", None, "°C (MOD11A2: factor 0.02, offset −273.15)",
     "Define el rango de operación diurna de equipos y trajes y el estrés térmico de materiales."),
    ("amplitud termica", "amplitud_c", None, "Δ°C (se calcula LST día − noche si falta)",
     "Los ciclos térmicos fuertes fatigan rocas, sellos y estructuras; es un proxy de inercia térmica baja."),
    ("irradiancia solar", "irradiancia_wm2", None, "W/m² (NASA POWER ALLSKY_SFC_SW_DWN en kWh/m²/día → ×41.667)",
     "Dimensiona paneles solares y balance térmico; en Marte llega ~43 % de la irradiancia terrestre."),
    ("opacidad atmosferica", "aod", None, "adimensional (AOD MODIS 550 nm)",
     "La carga de polvo reduce la energía solar y ensucia óptica; la AOD terrestre a 550 nm solo aproxima τ a 880 nm."),
    ("duracion del dia", "sol_h", None, "h",
     "El sol marciano (24.66 h) marca ritmo circadiano de la tripulación y ciclos de energía."),
    ("presion atmosferica", "presion_pa", None, "Pa",
     "La presión ~6 mbar decide estabilidad del agua líquida, paracaídas, convección y transporte de polvo."),
    ("densidad atmosferica", "densidad_atm", None, "kg/m³",
     "La atmósfera tenue limita el frenado aerodinámico en descenso y el enfriamiento convectivo."),
    ("velocidad del viento", "viento_ms", 10.0, "m/s (NASA POWER WS10M / WS2M)",
     "El viento mueve polvo, forma dunas y carga estructuras; en un análogo permite probar mitigación de polvo."),
    ("perdida de potencia", "polvo_pct_dia", 0.5, "%/día",
     "La acumulación de polvo en paneles reduce la energía disponible día con día (dato de los rovers MER)."),
    ("tormentas de polvo", "tormentas_polvo", None, "texto (local / regional) o MERRA-2",
     "Las tormentas de polvo condicionan la operación solar y la visibilidad; los análogos con eventos de polvo permiten ensayos."),
    ("formas y procesos eolicos", "formas_eolicas", None, "texto (dunas, ripples, yardangs, remolinos)",
     "Dunas, ripples y yardangs son geoformas marcianas comunes; sirven para entrenar interpretación de imágenes HiRISE y pruebas de movilidad."),
    ("contenido de agua en regolito", "agua_regolito_wt", 10.0, "wt% (SMAP es m³/m³ → ×100 aprox., volumétrico ≠ gravimétrico)",
     "El agua en el regolito es el recurso ISRU más importante (agua, O₂, propelente)."),
    ("profundidad del hielo", "hielo_prof_m", None, "m",
     "Hielo a <1 m es extraíble con excavación ligera; define sitios de base en Marte (mapa SWIM)."),
    ("fraccion volumetrica de hielo", "hielo_pct", None, "%vol",
     "La pureza del hielo subterráneo dimensiona la planta de extracción de agua."),
    ("tipo de recurso hidrico", "recurso_hidrico", None, "texto",
     "Clasifica el recurso de agua disponible (hielo enterrado, glaciar cubierto, minerales hidratados, permafrost)."),
    ("actividad de agua", "aw", 0.2, "adimensional",
     "aw < 0.6 deja el sitio fuera de 'región especial' (protección planetaria) y limita la vida microbiana."),
    ("ph de agua", "ph", 2.0, "pH",
     "El pH de salmueras y suelos condiciona habitabilidad, corrosión y química de extracción."),
    ("extremofilos", "extremofilos", None, "texto",
     "Comunidades extremófilas (halófilos, endolitos, acidófilos) son modelos de vida posible en Marte y de detección de biofirmas."),
    ("biomasa microbiana", "biomasa_cel_g", None, "células/g",
     "Una biomasa baja (tipo Atacama) reproduce el reto de detectar vida escasa con instrumentos de misión."),
    ("biofirmas", "biofirmas", None, "texto",
     "Biofirmas preservadas (lípidos, estromatolitos, δ¹³C, minerales biogénicos) permiten validar instrumentos de búsqueda de vida."),
    ("proteccion planetaria", "cospar", None, "texto",
     "La categoría COSPAR fija las reglas de limpieza y de búsqueda de vida en la misión."),
    ("tipo de cavidad", "cavidad", None, "texto",
     "Tubos de lava y claraboyas dan blindaje natural contra radiación y micrometeoritos para hábitats."),
    ("diametro interno", "tubo_diam_m", None, "m",
     "El diámetro del tubo determina si caben módulos de hábitat y vehículos."),
    ("espesor de techo", "techo_m", None, "m",
     "Unos metros de roca de techo bastan para blindar la radiación cósmica; también fijan la estabilidad."),
    ("dosis efectiva", "dosis_msv", None, "mSv/año",
     "La dosis de radiación es el principal límite de salud de la tripulación; en la Tierra no se reproduce."),
    ("estabilidad estructural", "estabilidad_techo", None, "texto",
     "Un techo estable (arco basáltico) es requisito para ocupar una cavidad."),
    ("retardo de comunicacion", "retardo_min", None, "min",
     "El retardo de 3–22 min obliga a operaciones autónomas; se simula en cualquier análogo."),
    ("error de navegacion", "error_nav_m", None, "m",
     "La precisión de navegación define el tamaño de la elipse y la separación entre activos."),
    ("potencia electrica", "potencia_kwe", None, "kWe",
     "Dimensiona la planta de energía del hábitat/spaceport."),
    ("resistencia a compresion", "compresion_mpa", None, "MPa (ensayo de laboratorio)",
     "El material local sinterizado o con azufre debe alcanzar 10–70 MPa para construir con ISRU."),
]


def conocimiento(nombre: str):
    n = norm(nombre)
    for clave, col, tol, conv, por_que in CONOCIMIENTO:
        if clave in n:
            return col, tol, conv, por_que
    return slug(nombre), None, "", ""


# ---------------------------------------------------------------------------
# Catálogo de variables
# ---------------------------------------------------------------------------
@dataclass
class Variable:
    id: str
    seccion: str
    nombre: str
    unidad_si: str
    unidad_en: str
    min_si: float | None
    max_si: float | None
    clases: list
    relevancia: float | None
    fidelidad: float | None
    fuente_marte: str
    fuente_analogo: str
    observacion: str
    tipo: str            # numerica | categorica | limite | mision | norma
    columna: str = ""
    tolerancia: float | None = None
    factor: float = 1.0
    offset: float = 0.0
    peso_indice: float | None = None
    conversion: str = ""
    por_que: str = ""


def leer_catalogo(path: Path):
    texto = leer_texto(path)
    filas = list(csv.reader(io.StringIO(texto)))
    variables, notas, seccion = [], [], ""
    encabezado_visto = False
    for f in filas:
        f = [reparar(c) for c in f] + [""] * 11
        c0 = f[0].strip()
        if not c0:
            continue
        if c0.lower().startswith("variable ["):
            encabezado_visto = True
            continue
        if not encabezado_visto:
            continue
        if c0.startswith("•"):
            notas.append(c0.lstrip("• ").strip())
            continue
        if all(not x.strip() for x in f[1:11]):
            if c0.lower() != "notas y supuestos":
                seccion = re.sub(r"\s*\(.*\)$", "", c0)
            continue
        nombre_full = c0
        m = re.search(r"\[([^\]|]*)\|([^\]]*)\]\s*$", nombre_full)
        u_si, u_en = (m.group(1).strip(), m.group(2).strip()) if m else ("", "")
        nombre = re.sub(r"\s*\[[^\]]*\]\s*$", "", nombre_full)
        if "—" in nombre:
            nombre = nombre.split("—", 1)[1].strip()
        clases = [c.strip() for c in f[5].split(";")] if f[5].strip() and f[5].strip() != "N/A" else []
        rel, fid = num(f[6]), num(f[7])
        fa = f[9]
        if fid is None:
            tipo = "norma"
        elif fid <= 1:
            tipo = "limite"
        elif re.search(r"no depende del sitio|se simula", norm(fa)):
            tipo = "mision"
        elif clases:
            tipo = "categorica"
        else:
            tipo = "numerica"
        col, tol, conv, por_que = conocimiento(nombre)
        v = Variable(
            id=slug(nombre), seccion=seccion, nombre=nombre, unidad_si=u_si, unidad_en=u_en,
            min_si=num(f[1]), max_si=num(f[2]), clases=clases, relevancia=rel, fidelidad=fid,
            fuente_marte=f[8], fuente_analogo=fa, observacion=f[10], tipo=tipo,
            columna=col, tolerancia=tol, conversion=conv, por_que=por_que,
        )
        m30 = re.search(r"pesa\s*(\d+)\s*%", f[10])
        if m30:
            v.peso_indice = float(m30.group(1)) / 100
        variables.append(v)
    return variables, notas


def crear_mapeo(variables, destino: Path):
    with destino.open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["id_variable", "variable", "seccion", "tipo", "columna_capa", "factor", "offset",
                    "tolerancia", "peso_indice", "unidad_destino_SI", "nota_conversion"])
        for v in variables:
            w.writerow([v.id, v.nombre, v.seccion, v.tipo, v.columna, 1, 0,
                        "" if v.tolerancia is None else v.tolerancia,
                        "" if v.peso_indice is None else v.peso_indice, v.unidad_si, v.conversion])
    print(f"✓ Plantilla de mapeo creada: {destino}\n  Edita la columna 'columna_capa' con el nombre real del atributo en tu capa QGIS.")


def aplicar_mapeo(variables, path: Path):
    filas = list(csv.DictReader(io.StringIO(leer_texto(path))))
    por_id = {r["id_variable"].strip(): r for r in filas}
    for v in variables:
        r = por_id.get(v.id)
        if not r:
            continue
        v.columna = (r.get("columna_capa") or "").strip()
        v.factor = num(r.get("factor")) or 1.0
        v.offset = num(r.get("offset")) or 0.0
        t = num(r.get("tolerancia"))
        v.tolerancia = t if t is not None else v.tolerancia
        p = num(r.get("peso_indice"))
        v.peso_indice = p if p is not None else v.peso_indice


# ---------------------------------------------------------------------------
# Lectura de la capa de zonas (GeoJSON / GPKG nativo, o geopandas)
# ---------------------------------------------------------------------------
def _wkb_polys(b: bytes, o: int = 0):
    """Devuelve (lista de polígonos [anillos [(x,y)]], offset)."""
    bo = "<" if b[o] == 1 else ">"
    t = struct.unpack_from(bo + "I", b, o + 1)[0]
    o += 5
    if t & 0x20000000:          # EWKB con SRID
        o += 4
    base = t & 0xFFFF
    dims = 2
    if base > 3000:
        base -= 3000; dims = 4
    elif base > 2000:
        base -= 2000; dims = 3
    elif base > 1000:
        base -= 1000; dims = 3
    if t & 0x80000000 or t & 0x40000000:
        dims = 2 + bool(t & 0x80000000) + bool(t & 0x40000000)

    def ring(o):
        n = struct.unpack_from(bo + "I", b, o)[0]; o += 4
        pts = []
        for _ in range(n):
            x, y = struct.unpack_from(bo + "dd", b, o)
            pts.append((x, y)); o += 8 * dims
        return pts, o

    if base == 3:
        nr = struct.unpack_from(bo + "I", b, o)[0]; o += 4
        rings = []
        for _ in range(nr):
            r, o = ring(o); rings.append(r)
        return [rings], o
    if base == 6:
        n = struct.unpack_from(bo + "I", b, o)[0]; o += 4
        polys = []
        for _ in range(n):
            p, o = _wkb_polys(b, o); polys += p
        return polys, o
    if base == 1:
        x, y = struct.unpack_from(bo + "dd", b, o)
        return [[[(x, y)]]], o + 8 * dims
    raise ValueError(f"Tipo de geometría WKB no soportado: {base}")


def _gpkg_geom(blob: bytes):
    if blob is None:
        return []
    if blob[:2] != b"GP":
        return _wkb_polys(blob)[0]
    flags = blob[3]
    env = (flags >> 1) & 0x07
    env_len = {0: 0, 1: 32, 2: 48, 3: 48, 4: 64}.get(env, 0)
    return _wkb_polys(blob, 8 + env_len)[0]


def leer_gpkg(path: Path, capa: str | None):
    con = sqlite3.connect(str(path))
    cur = con.cursor()
    cur.execute("SELECT table_name, column_name, srs_id FROM gpkg_geometry_columns")
    capas = cur.fetchall()
    if not capas:
        raise SystemExit("El GeoPackage no tiene capas vectoriales.")
    tabla, gcol, srs = next((c for c in capas if c[0] == capa), capas[0])
    if srs not in (4326, 0, -1):
        print(f"⚠ La capa '{tabla}' está en EPSG:{srs}. Sin geopandas no se reproyecta: "
              "el centroide y el área se omiten. Exporta en EPSG:4326 desde QGIS o instala geopandas.")
    cur.execute(f'SELECT * FROM "{tabla}"')
    cols = [d[0] for d in cur.description]
    zonas = []
    for fila in cur.fetchall():
        props = dict(zip(cols, fila))
        geom = props.pop(gcol)
        polys = _gpkg_geom(geom) if srs in (4326, 0, -1) else []
        zonas.append({"props": props, "polys": polys})
    return zonas


def leer_geojson(path: Path):
    gj = json.loads(leer_texto(path))
    zonas = []
    for ft in gj.get("features", []):
        g = ft.get("geometry") or {}
        t, c = g.get("type"), g.get("coordinates")
        if t == "Polygon":
            polys = [[[tuple(p[:2]) for p in r] for r in c]]
        elif t == "MultiPolygon":
            polys = [[[tuple(p[:2]) for p in r] for r in poly] for poly in c]
        elif t == "Point":
            polys = [[[tuple(c[:2])]]]
        else:
            polys = []
        zonas.append({"props": ft.get("properties") or {}, "polys": polys})
    return zonas


def leer_csv_zonas(path: Path):
    filas = list(csv.DictReader(io.StringIO(leer_texto(path))))
    zonas = []
    for r in filas:
        lat = num(r.get("lat") or r.get("latitud") or r.get("y"))
        lon = num(r.get("lon") or r.get("longitud") or r.get("x"))
        polys = [[[(lon, lat)]]] if lat is not None and lon is not None else []
        zonas.append({"props": r, "polys": polys})
    return zonas


def leer_zonas(path: Path, capa: str | None = None):
    suf = path.suffix.lower()
    try:
        import geopandas as gpd  # type: ignore
        if suf != ".csv":
            gdf = gpd.read_file(path, layer=capa) if capa else gpd.read_file(path)
            if gdf.crs is not None and gdf.crs.to_epsg() != 4326:
                gdf = gdf.to_crs(4326)
            gj = json.loads(gdf.to_json())
            tmp = Path(".__zonas_tmp.geojson")
            tmp.write_text(json.dumps(gj), encoding="utf-8")
            z = leer_geojson(tmp); tmp.unlink()
            return z
    except ImportError:
        pass
    if suf in (".geojson", ".json"):
        return leer_geojson(path)
    if suf == ".gpkg":
        return leer_gpkg(path, capa)
    if suf == ".csv":
        return leer_csv_zonas(path)
    raise SystemExit(f"Formato {suf} requiere geopandas. Exporta la capa a GeoJSON o GeoPackage desde QGIS.")


def geometria_resumen(polys):
    """Centroide, área (km²), extensión máxima (km) y bbox, con proyección local."""
    pts = [p for poly in polys for r in poly for p in r]
    if not pts:
        return None
    lat0 = sum(p[1] for p in pts) / len(pts)
    kx = 111.320 * math.cos(math.radians(lat0)); ky = 110.574
    A = cx = cy = 0.0
    for poly in polys:
        for i, r in enumerate(poly):
            if len(r) < 3:
                continue
            s = 0.0; sx = sy = 0.0
            for (x1, y1), (x2, y2) in zip(r, r[1:] + r[:1]):
                X1, Y1, X2, Y2 = x1 * kx, y1 * ky, x2 * kx, y2 * ky
                c = X1 * Y2 - X2 * Y1
                s += c; sx += (X1 + X2) * c; sy += (Y1 + Y2) * c
            s /= 2
            sign = 1 if i == 0 else -1
            a = abs(s) * sign
            if s != 0:
                cx += sx / (6 * s) * a; cy += sy / (6 * s) * a
            A += a
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    if A > 0:
        clon, clat = cx / A / kx, cy / A / ky
    else:
        clon, clat = sum(xs) / len(xs), sum(ys) / len(ys)
    ext = math.hypot((max(xs) - min(xs)) * kx, (max(ys) - min(ys)) * ky)
    return {"lat": clat, "lon": clon, "area_km2": A, "extension_km": ext,
            "bbox": [min(xs), min(ys), max(xs), max(ys)]}


# ---------------------------------------------------------------------------
# Evaluación de cada variable
# ---------------------------------------------------------------------------
def buscar_col(props: dict, col: str):
    if not col:
        return None
    if col in props:
        return props[col]
    nc = norm(col)
    for k, v in props.items():
        if norm(k) == nc:
            return v
    return None


def puntuar_rango(v, lo, hi, tol):
    if lo is None and hi is None:
        return None, "sin_rango", None
    lo = -math.inf if lo is None else lo
    hi = math.inf if hi is None else hi
    if lo <= v <= hi:
        return 1.0, "dentro", 0.0
    d = (lo - v) if v < lo else (v - hi)
    ancho = (hi - lo) if math.isfinite(hi - lo) else 0
    tol = tol or max(ancho, abs(lo) if math.isfinite(lo) else 0, abs(hi) if math.isfinite(hi) else 0, 1e-6)
    s = max(0.0, 1 - d / tol)
    return s, ("debajo" if v < lo else "encima"), d


def palabras_clave(clase: str):
    base = norm(re.sub(r"\(.*?\)", "", clase))
    toks = [t for t in re.split(r"[ /,]+", base) if len(t) >= 5 and t not in
            ("latitudes", "medias", "cubiertos", "escombros", "busqueda", "regiones")]
    return base, toks


def evaluar_categoria(valor: str, clases):
    nv = norm(valor)
    hits = []
    for c in clases:
        base, toks = palabras_clave(c)
        if base and base in nv:
            hits.append(c)
        elif any(t in nv for t in toks):
            hits.append(c)
    return (1.0 if hits else 0.0), hits


def ESTADO(score):
    if score is None:
        return "sin dato"
    if score >= 0.999:
        return "dentro"
    if score >= 0.5:
        return "cercano"
    return "fuera"


def evaluar_zona(z, variables, geo):
    props = z["props"]
    res = []
    for v in variables:
        r = {"id": v.id, "seccion": v.seccion, "variable": v.nombre, "tipo": v.tipo,
             "unidad": v.unidad_si, "min_marte": v.min_si, "max_marte": v.max_si,
             "clases_marte": v.clases, "relevancia": v.relevancia, "fidelidad": v.fidelidad,
             "fuente_marte": v.fuente_marte, "fuente_analogo": v.fuente_analogo,
             "observacion": v.observacion, "por_que": v.por_que, "columna": v.columna,
             "valor": None, "score": None, "estado": "sin dato", "detalle": "", "coincidencias": []}
        bruto = buscar_col(props, v.columna)
        # valores derivados automáticamente
        if bruto in (None, "") and v.id.startswith("amplitud_termica"):
            d = buscar_col(props, next((x.columna for x in variables if x.id.startswith("temperatura_maxima")), ""))
            n = buscar_col(props, next((x.columna for x in variables if x.id.startswith("temperatura_minima")), ""))
            vd = [x for x in variables if x.id.startswith("temperatura_maxima")][0]
            vn = [x for x in variables if x.id.startswith("temperatura_minima")][0]
            if num(d) is not None and num(n) is not None:
                bruto = (num(d) * vd.factor + vd.offset) - (num(n) * vn.factor + vn.offset)
                r["detalle"] = "Calculada: LST día − LST noche. "
        if bruto in (None, "") and v.id.startswith("semieje_mayor") and geo:
            bruto = geo["extension_km"]; r["detalle"] = "Calculada: extensión máxima del polígono. "
        if bruto in (None, "") and v.id.startswith("latitud_del_sitio") and geo:
            bruto = geo["lat"]; r["detalle"] = "Calculada: latitud del centroide. "

        if v.tipo == "categorica":
            if bruto not in (None, ""):
                s, hits = evaluar_categoria(str(bruto), v.clases)
                r.update(valor=str(bruto), score=s, coincidencias=hits, estado="dentro" if s else "fuera")
        elif v.tipo in ("numerica", "limite", "mision"):
            x = num(bruto)
            if x is not None:
                if not r["detalle"].startswith("Calculada: LST"):
                    x = x * v.factor + v.offset
                r["valor"] = x
                if v.tipo == "numerica":
                    if v.id.startswith("semieje_mayor"):
                        s = 1.0 if x >= (v.min_si or 0) else x / (v.min_si or 1)
                        r["score"], r["estado"] = s, ESTADO(s)
                        r["detalle"] += f"Cabe una elipse de hasta {fmt(min(x, v.max_si or x))} km."
                    else:
                        s, est, d = puntuar_rango(x, v.min_si, v.max_si, v.tolerancia)
                        r["score"], r["estado"] = s, ESTADO(s)
                        if d:
                            r["detalle"] += f"{'Por debajo' if est == 'debajo' else 'Por encima'} de la ventana marciana por {fmt(d)} {v.unidad_si}."
        else:  # norma
            if bruto not in (None, ""):
                r["valor"] = str(bruto)
        res.append(r)
    return res


def peso(r):
    return (r["relevancia"] or 0) * (r["fidelidad"] or 0)


def indices(res, variables):
    evaluables = [r for r in res if r["tipo"] in ("numerica", "categorica")]
    con_dato = [r for r in evaluables if r["score"] is not None]
    w_total = sum(peso(r) for r in evaluables) or 1
    w_dato = sum(peso(r) for r in con_dato)
    fijos = {v.id: v.peso_indice for v in variables if v.peso_indice}
    fijo_presentes = {k: p for k, p in fijos.items() if any(r["id"] == k and r["score"] is not None for r in res)}
    resto = [r for r in con_dato if r["id"] not in fijos]
    w_resto = sum(peso(r) for r in resto)
    s_resto = sum(peso(r) * r["score"] for r in resto) / w_resto if w_resto else None
    peso_resto = max(0.0, 1 - sum(fijos.values()))
    num_, den = 0.0, 0.0
    for k, p in fijo_presentes.items():
        num_ += p * next(r["score"] for r in res if r["id"] == k); den += p
    if s_resto is not None:
        num_ += peso_resto * s_resto; den += peso_resto
    indice = num_ / den if den else None
    secciones = {}
    for r in con_dato:
        s = secciones.setdefault(r["seccion"], [0.0, 0.0, 0])
        s[0] += peso(r) * r["score"]; s[1] += peso(r); s[2] += 1
    sec = {k: (a / b if b else None, n) for k, (a, b, n) in secciones.items()}
    return {"indice": indice, "cobertura": w_dato / w_total, "n_evaluadas": len(con_dato),
            "n_evaluables": len(evaluables), "secciones": sec,
            "indice_ponderado_simple": (sum(peso(r) * r["score"] for r in con_dato) / w_dato) if w_dato else None}


def clase_analogia(i, cob):
    if i is None:
        return "Sin datos suficientes"
    c = "Análogo alto" if i >= 0.8 else "Análogo moderado" if i >= 0.6 else "Análogo parcial" if i >= 0.4 else "Análogo débil"
    if cob < 0.4:
        c += " (confianza baja: pocas variables medidas)"
    return c


USOS = [
    ({"Geología", "Topografía"}, "Pruebas de aterrizaje, movilidad de rovers y cimentación de hábitats/plataformas"),
    ({"Clima", "Temperatura e iluminación"}, "Ensayos de operación en hiperaridez, gestión térmica y energía solar"),
    ({"Atmósfera, polvo y viento"}, "Pruebas de mitigación de polvo en paneles, sellos y óptica"),
    ({"Agua y recursos"}, "Prospección y extracción de agua (ISRU)"),
    ({"Biología y habitabilidad"}, "Astrobiología: detección de biofirmas y calibración de instrumentos de búsqueda de vida"),
    ({"Subsuelo y cuevas"}, "Hábitats en cavidades y blindaje natural contra radiación"),
]


# ---------------------------------------------------------------------------
# Criterios de spaceport terrestre (docx "Variables por definir")
# ---------------------------------------------------------------------------
def leer_spaceport(path: Path | None):
    if not path or not path.exists():
        return []
    return list(csv.DictReader(io.StringIO(leer_texto(path))))


def evaluar_spaceport(props, criterios, geo):
    out = []
    for c in criterios:
        col = (c.get("columna_capa") or "").strip()
        bruto = buscar_col(props, col)
        if bruto in (None, "") and col == "lat_abs" and geo:
            bruto = abs(geo["lat"])
        x = num(bruto)
        lo, hi = num(c.get("min")), num(c.get("max"))
        r = {"parametro": c["parametro"], "unidad": c.get("unidad", ""), "min": lo, "max": hi,
             "ideal": c.get("ideal", ""), "utilidad": c.get("utilidad", ""), "alternativa": c.get("alternativa", ""),
             "valor": x if x is not None else (bruto if bruto not in (None, "") else None),
             "score": None, "estado": "sin dato"}
        if x is not None and (lo is not None or hi is not None):
            s, est, d = puntuar_rango(x, lo, hi, num(c.get("tolerancia")))
            r["score"], r["estado"] = s, ESTADO(s)
        elif lo is None and hi is None:
            r["estado"] = "cualitativo"
        out.append(r)
    con = [r for r in out if r["score"] is not None]
    idx = sum(r["score"] for r in con) / len(con) if con else None
    return out, idx


# ---------------------------------------------------------------------------
# Redacción de la ficha
# ---------------------------------------------------------------------------
def rango_txt(r):
    if r["tipo"] == "categorica":
        return "; ".join(r["clases_marte"])
    lo, hi = r["min_marte"], r["max_marte"]
    if lo is None and hi is None:
        return "—"
    if lo == hi:
        return f"{fmt(lo)} {r['unidad']}"
    return f"{fmt(lo)} a {fmt(hi)} {r['unidad']}"


def valor_txt(r):
    if r["valor"] is None:
        return "—"
    if isinstance(r["valor"], str):
        return r["valor"]
    return f"{fmt(r['valor'])} {r['unidad']}"


def lc(s):
    return s[:1].lower() + s[1:] if s[:2].isupper() is False else s


def narrativa(nombre, res, idx, geo):
    eval_ = [r for r in res if r["score"] is not None and r["tipo"] in ("numerica", "categorica")]
    fuertes = sorted([r for r in eval_ if r["score"] >= 0.999], key=lambda r: -peso(r))
    cercanas = sorted([r for r in eval_ if 0.5 <= r["score"] < 0.999], key=lambda r: -peso(r))
    debiles = sorted([r for r in eval_ if r["score"] < 0.5], key=lambda r: -peso(r))
    p = []
    i, cob = idx["indice"], idx["cobertura"]
    if i is None:
        p.append(f"{nombre} aún no tiene variables medidas suficientes para calcular el índice de analogía.")
    else:
        p.append(f"{nombre} obtiene un índice de analogía con Marte de {i:.2f} (escala 0–1), "
                 f"calculado con {idx['n_evaluadas']} de {idx['n_evaluables']} variables evaluables "
                 f"({cob*100:.0f} % del peso total). Se clasifica como **{clase_analogia(i, cob)}**.")
    if fuertes:
        items = []
        for r in fuertes[:6]:
            if r["id"].startswith("semieje_mayor"):
                items.append(f"extensión de {fmt(r['valor'])} km, suficiente para trazar una elipse de aterrizaje de {rango_txt(r)}")
            elif r["tipo"] == "categorica":
                items.append(f"{lc(r['variable'])} ({', '.join(r['coincidencias'])}, igual que en Marte)")
            else:
                items.append(f"{lc(r['variable'])} de {valor_txt(r)}, dentro de la ventana marciana ({rango_txt(r)})")
        p.append("Es análogo principalmente por: " + "; ".join(items) + ".")
        razones = [r["por_que"] for r in fuertes[:3] if r["por_que"]]
        if razones:
            p.append("Esto importa porque: " + " ".join(razones))
    if cercanas:
        p.append("Se acerca, sin igualarlas, en: " + "; ".join(
            f"{lc(r['variable'])} ({valor_txt(r)} frente a {rango_txt(r)})" for r in cercanas[:4]) + ".")
    if debiles:
        p.append("Difiere de Marte en: " + "; ".join(
            f"{lc(r['variable'])} ({valor_txt(r)} frente a {rango_txt(r)})" for r in debiles[:4]) +
            ". Estas variables deben declararse como límites del análogo o compensarse en laboratorio.")
    usos = []
    for secs, uso in USOS:
        vals = [idx["secciones"][s][0] for s in secs if s in idx["secciones"] and idx["secciones"][s][0] is not None]
        n = sum(idx["secciones"][s][1] for s in secs if s in idx["secciones"])
        if vals and n >= 2 and sum(vals) / len(vals) >= 0.7:
            usos.append(uso)
    return p, usos


def construir_ficha(z, variables, criterios, id_col, nombre_col):
    props = z["props"]
    geo = geometria_resumen(z["polys"])
    zid = buscar_col(props, id_col) if id_col else None
    if zid in (None, ""):
        zid = props.get("fid") or props.get("id_zona") or ""
    zid = str(zid)
    nombre = buscar_col(props, nombre_col) or props.get("nombre") or props.get("name") or f"Zona {zid}"
    nombre = str(nombre)
    if not zid:
        zid = slug(nombre)
    zid = re.sub(r"[^A-Za-z0-9_\-]+", "_", zid)
    res = evaluar_zona(z, variables, geo)
    idx = indices(res, variables)
    parrafos, usos = narrativa(nombre, res, idx, geo)
    sp, sp_idx = evaluar_spaceport(props, criterios, geo)
    usadas = {v.columna for v in variables} | {c.get("columna_capa") for c in criterios} | {id_col, nombre_col}
    extra = {k: v for k, v in props.items() if k not in usadas and not isinstance(v, (bytes, bytearray))}
    return {
        "id": zid, "nombre": nombre, "geo": geo, "indice": idx["indice"], "cobertura": idx["cobertura"],
        "clase": clase_analogia(idx["indice"], idx["cobertura"]), "secciones": idx["secciones"],
        "n_evaluadas": idx["n_evaluadas"], "n_evaluables": idx["n_evaluables"],
        "narrativa": parrafos, "usos_recomendados": usos, "variables": res,
        "spaceport": sp, "indice_spaceport": sp_idx, "atributos_extra": extra,
        "polys": z["polys"],
    }


# ---------------------------------------------------------------------------
# Salidas: Markdown, HTML, JSON, CSV, GeoJSON, KML
# ---------------------------------------------------------------------------
ICON = {"dentro": "✔", "cercano": "≈", "fuera": "✘", "sin dato": "—", "cualitativo": "◦"}


def ficha_md(f, notas):
    L = [f"# Ficha científica — {f['nombre']}", ""]
    g = f["geo"]
    L += [f"**ID:** {f['id']}  ", f"**Índice de analogía con Marte:** {fmt(f['indice'])} — {f['clase']}  ",
          f"**Cobertura de datos:** {f['cobertura']*100:.0f} % ({f['n_evaluadas']}/{f['n_evaluables']} variables)  "]
    if g:
        L.append(f"**Centroide:** {g['lat']:.5f}, {g['lon']:.5f} · **Área:** {fmt(g['area_km2'])} km² · **Extensión:** {fmt(g['extension_km'])} km  ")
    L += ["", "## ¿Por qué es análogo de Marte?", ""] + [x + "\n" for x in f["narrativa"]]
    if f["usos_recomendados"]:
        L += ["**Usos recomendados del análogo:**", ""] + [f"- {u}" for u in f["usos_recomendados"]] + [""]
    L += ["## Puntuación por sección", "", "| Sección | Puntuación | Variables con dato |", "|---|---|---|"]
    for s, (sc, n) in f["secciones"].items():
        L.append(f"| {s} | {fmt(sc)} | {n} |")
    L.append("")
    sec = None
    L += ["## Propiedades del sitio frente a la ventana marciana", ""]
    for r in [r for r in f["variables"] if r["tipo"] in ("numerica", "categorica")]:
        if r["seccion"] != sec:
            sec = r["seccion"]
            L += ["", f"### {sec}", "", "| Variable | Valor del sitio | Ventana Marte | Estado | Punt. | R | F | Fuente Marte | Fuente análogo |",
                  "|---|---|---|---|---|---|---|---|---|"]
        L.append(f"| {r['variable']} | {valor_txt(r)} | {rango_txt(r)} | {ICON[r['estado']]} {r['estado']} | {fmt(r['score'])} | "
                 f"{fmt(r['relevancia'])} | {fmt(r['fidelidad'])} | {r['fuente_marte']} | {r['fuente_analogo']} |")
    L += ["", "## Lo que la Tierra no puede imitar (límites de la analogía)", "",
          "| Variable | Valor en Marte | Valor del sitio | Observación |", "|---|---|---|---|"]
    for r in [r for r in f["variables"] if r["tipo"] == "limite"]:
        L.append(f"| {r['variable']} | {rango_txt(r)} | {valor_txt(r)} | {r['observacion'] or r['fuente_analogo']} |")
    falt = [r for r in f["variables"] if r["tipo"] in ("numerica", "categorica") and r["score"] is None]
    if falt:
        L += ["", "## Pendiente de medir (campo o laboratorio)", ""]
        L += [f"- **{r['variable']}** ({rango_txt(r)}) → {r['fuente_analogo']}" for r in falt]
    L += ["", "## Parámetros de misión (no dependen del sitio)", ""]
    L += [f"- {r['variable']}: {rango_txt(r)} — {r['fuente_marte']}" for r in f["variables"] if r["tipo"] in ("mision", "norma")]
    if f["spaceport"]:
        L += ["", f"## Aptitud como spaceport terrestre (índice {fmt(f['indice_spaceport'])})", "",
              "| Parámetro | Valor | Criterio | Estado | Utilidad |", "|---|---|---|---|---|"]
        for r in f["spaceport"]:
            crit = r["ideal"] or f"{fmt(r['min'])} a {fmt(r['max'])} {r['unidad']}"
            val = fmt(r["valor"]) if not isinstance(r["valor"], str) else r["valor"]
            L.append(f"| {r['parametro']} | {val} | {crit} | {ICON[r['estado']]} {r['estado']} | {r['utilidad']} |")
    if f["atributos_extra"]:
        L += ["", "## Otros atributos de la capa", ""] + [f"- {k}: {v}" for k, v in f["atributos_extra"].items()]
    L += ["", "## Notas metodológicas", ""] + [f"- {n}" for n in notas] + [
        "- Puntuación por variable: 1 si el valor cae dentro de la ventana marciana; fuera, decae linealmente con la distancia (tolerancia = ancho de la ventana o la definida en mapeo_columnas.csv). Categóricas: 1 si la clase observada coincide con una clase marciana.",
        "- Índice de analogía: precipitación 30 % + NDVI 30 % + 40 % media de las demás variables ponderada por relevancia × fidelidad. Las variables de fidelidad 1 no puntúan: se reportan como límites."]
    return "\n".join(L)


CSS = """
:root{--ink:#1d1d1f;--mut:#6b6b70;--line:#e4e1dc;--bg:#faf8f5;--mars:#b5462b;--ok:#2f7d4f;--mid:#b07a12;--bad:#a8322a}
*{box-sizing:border-box}body{margin:0;font:15px/1.55 -apple-system,"Segoe UI",Inter,Roboto,sans-serif;color:var(--ink);background:var(--bg)}
.wrap{max-width:1100px;margin:0 auto;padding:40px 28px 80px}
header{border-bottom:3px solid var(--mars);padding-bottom:18px;margin-bottom:26px}
.kick{text-transform:uppercase;letter-spacing:.14em;font-size:12px;color:var(--mars);font-weight:700}
h1{font-size:34px;margin:6px 0 10px;line-height:1.15}h2{font-size:20px;margin:36px 0 12px;border-bottom:1px solid var(--line);padding-bottom:6px}
h3{font-size:15px;margin:22px 0 8px;color:var(--mars);text-transform:uppercase;letter-spacing:.06em}
.meta{display:flex;flex-wrap:wrap;gap:10px 26px;color:var(--mut);font-size:14px}
.kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:12px;margin:22px 0}
.kpi{background:#fff;border:1px solid var(--line);border-radius:10px;padding:14px 16px}
.kpi b{display:block;font-size:26px}.kpi span{font-size:12px;color:var(--mut);text-transform:uppercase;letter-spacing:.06em}
.bar{height:6px;background:#eee;border-radius:3px;margin-top:8px;overflow:hidden}.bar i{display:block;height:100%;background:var(--mars)}
.lead p{margin:0 0 10px}.usos li{margin:3px 0}
table{width:100%;border-collapse:collapse;font-size:13px;background:#fff;border:1px solid var(--line)}
th,td{padding:7px 9px;border-bottom:1px solid var(--line);vertical-align:top;text-align:left}
th{background:#f3f0ea;font-weight:600;font-size:12px}td.n{white-space:nowrap}
.st{font-weight:600;white-space:nowrap}.dentro{color:var(--ok)}.cercano{color:var(--mid)}.fuera{color:var(--bad)}.sin{color:var(--mut)}
.src{color:var(--mut);font-size:12px}.why{color:var(--mut);font-size:12px;margin-top:3px}
.note{font-size:13px;color:var(--mut)}.pill{display:inline-block;padding:2px 9px;border-radius:20px;background:#f2e4df;color:var(--mars);font-weight:600;font-size:13px}
footer{margin-top:40px;font-size:12px;color:var(--mut)}
@media print{body{background:#fff}.wrap{padding:0}h2{break-after:avoid}tr{break-inside:avoid}}
"""


def e(x):
    return html.escape("" if x is None else str(x))


FICHA_CSS = r"""
@import url('https://fonts.googleapis.com/css2?family=Space+Grotesk:wght@400;500;600;700&family=Inter:wght@400;500;600&family=JetBrains+Mono:wght@400;600&display=swap');
:root{
  --bg:#07080a;--bg2:#0d0f13;--card:#12151b;--card2:#171b22;--line:#232833;--line2:#2e3441;
  --ink:#f3f1ee;--ink2:#c9c6c1;--mut:#878b94;--mut2:#5c616b;
  --mars:#c1440e;--mars2:#b5462b;--burn:#e0782f;--glow:rgba(193,68,14,.35);
  --ok:#3ecf8e;--okbg:rgba(62,207,142,.12);--mid:#f2b134;--midbg:rgba(242,177,52,.12);
  --bad:#f0573f;--badbg:rgba(240,87,63,.12);--na:#7d828c;--nabg:rgba(125,130,140,.14);
  --r:16px;--disp:'Space Grotesk',system-ui,-apple-system,'Segoe UI',sans-serif;
  --body:'Inter',system-ui,-apple-system,'Segoe UI',Roboto,sans-serif;--mono:'JetBrains Mono',ui-monospace,Menlo,monospace;
}
*{box-sizing:border-box}html{scroll-behavior:smooth;scroll-padding-top:80px}
body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.6 var(--body);-webkit-font-smoothing:antialiased}
a{color:inherit}
.wrap{max-width:1180px;margin:0 auto;padding:0 28px}
/* ---------- barra superior ---------- */
.top{position:sticky;top:0;z-index:40;background:rgba(7,8,10,.72);backdrop-filter:blur(14px);-webkit-backdrop-filter:blur(14px);border-bottom:1px solid var(--line)}
.top .wrap{display:flex;align-items:center;gap:22px;height:58px}
.brand{font:600 12px/1 var(--mono);letter-spacing:.18em;text-transform:uppercase;color:var(--ink2);white-space:nowrap;display:flex;align-items:center;gap:10px}
.brand i{width:9px;height:9px;border-radius:50%;background:var(--mars);box-shadow:0 0 12px var(--mars)}
.nav{display:flex;gap:4px;margin-left:auto;overflow-x:auto;scrollbar-width:none}.nav::-webkit-scrollbar{display:none}
.nav a{font:500 13px var(--body);color:var(--mut);text-decoration:none;padding:7px 12px;border-radius:999px;white-space:nowrap;transition:.2s}
.nav a:hover,.nav a.on{color:var(--ink);background:var(--card2)}
/* ---------- hero ---------- */
.hero{position:relative;overflow:hidden;padding:84px 0 56px;border-bottom:1px solid var(--line)}
.hero:before{content:"";position:absolute;inset:-30% -10% auto auto;width:900px;height:900px;border-radius:50%;
  background:radial-gradient(circle at 35% 35%,#e0782f 0%,#c1440e 22%,#6e2109 48%,rgba(30,10,4,0) 70%);opacity:.38;filter:blur(8px);pointer-events:none}
.hero:after{content:"";position:absolute;inset:0;pointer-events:none;
  background-image:radial-gradient(1px 1px at 12% 22%,#fff8 50%,transparent 51%),radial-gradient(1px 1px at 72% 12%,#fff6 50%,transparent 51%),
  radial-gradient(1.4px 1.4px at 38% 70%,#fff5 50%,transparent 51%),radial-gradient(1px 1px at 88% 82%,#fff7 50%,transparent 51%),
  radial-gradient(1px 1px at 22% 88%,#fff4 50%,transparent 51%),radial-gradient(1.2px 1.2px at 56% 40%,#fff3 50%,transparent 51%);
  background-size:520px 520px}
.hero .wrap{position:relative;z-index:1}
.hgrid{display:grid;grid-template-columns:1.5fr 1fr;gap:48px;align-items:center}
.kick{font:600 12px var(--mono);letter-spacing:.22em;text-transform:uppercase;color:var(--burn);display:flex;align-items:center;gap:12px}
.kick:before{content:"";width:28px;height:1px;background:var(--burn)}
h1{font:700 clamp(36px,6vw,68px)/1.02 var(--disp);letter-spacing:-.025em;margin:18px 0 22px;text-wrap:balance}
.chips{display:flex;flex-wrap:wrap;gap:8px}
.chip{font:500 12.5px var(--mono);color:var(--ink2);background:rgba(255,255,255,.04);border:1px solid var(--line2);padding:6px 11px;border-radius:999px}
.cls{display:inline-flex;align-items:center;gap:9px;margin-top:22px;font:600 15px var(--disp);padding:10px 18px 10px 14px;border-radius:999px;
  background:linear-gradient(90deg,rgba(193,68,14,.28),rgba(193,68,14,.08));border:1px solid rgba(224,120,47,.45);color:#ffd9c4}
.cls i{width:8px;height:8px;border-radius:50%;background:var(--burn);box-shadow:0 0 10px var(--burn);animation:pulse 2.4s infinite}
@keyframes pulse{50%{opacity:.35}}
/* indicador circular */
.gauge{justify-self:end;position:relative;width:300px;max-width:100%;aspect-ratio:1}
.gauge svg{width:100%;height:100%;transform:rotate(-90deg)}
.gauge .trk{stroke:rgba(255,255,255,.07)}.gauge .val{stroke:url(#gg);stroke-linecap:round;transition:stroke-dashoffset 1.8s cubic-bezier(.2,.8,.2,1);filter:drop-shadow(0 0 10px var(--glow))}
.gauge .ctr{position:absolute;inset:0;display:flex;flex-direction:column;align-items:center;justify-content:center;text-align:center}
.gauge .big{font:700 74px/1 var(--disp);letter-spacing:-.04em}
.gauge .lab{font:600 11px var(--mono);letter-spacing:.2em;text-transform:uppercase;color:var(--mut);margin-top:10px}
.gauge .sub{font:500 12px var(--mono);color:var(--ink2);margin-top:4px}
/* KPI */
.kpis{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin-top:56px}
.kpi{background:linear-gradient(180deg,var(--card2),var(--card));border:1px solid var(--line);border-radius:var(--r);padding:22px 22px 20px;transition:.25s;position:relative;overflow:hidden}
.kpi:hover{border-color:rgba(224,120,47,.5);transform:translateY(-3px);box-shadow:0 18px 40px -20px var(--glow)}
.kpi .l{font:600 11px var(--mono);letter-spacing:.16em;text-transform:uppercase;color:var(--mut)}
.kpi .v{font:700 38px/1.1 var(--disp);letter-spacing:-.02em;margin-top:10px}.kpi .v small{font-size:16px;color:var(--mut);font-weight:500;margin-left:4px}
.kpi .d{font-size:12.5px;color:var(--mut);margin-top:4px}
.bar{height:6px;background:rgba(255,255,255,.07);border-radius:9px;overflow:hidden;margin-top:14px}
.bar>i{display:block;height:100%;width:0;border-radius:9px;background:linear-gradient(90deg,var(--mars),var(--burn));transition:width 1.4s cubic-bezier(.2,.8,.2,1)}
.bar.ok>i{background:var(--ok)}.bar.mid>i{background:var(--mid)}.bar.bad>i{background:var(--bad)}
/* conteo de estados */
.tally{display:flex;flex-wrap:wrap;gap:10px;margin-top:22px}
.tbtn{display:flex;align-items:center;gap:10px;font:500 13px var(--body);color:var(--ink2);background:var(--card);border:1px solid var(--line);padding:9px 14px;border-radius:12px;cursor:pointer;transition:.2s}
.tbtn b{font:700 18px var(--disp);color:var(--ink)}
.tbtn:hover{border-color:var(--line2);background:var(--card2)}
.tbtn.on{border-color:currentColor}
.tbtn[data-f=dentro]{--c:var(--ok)}.tbtn[data-f=cercano]{--c:var(--mid)}.tbtn[data-f=fuera]{--c:var(--bad)}.tbtn[data-f="sin dato"]{--c:var(--na)}
.tbtn .dot{width:10px;height:10px;border-radius:50%;background:var(--c);box-shadow:0 0 10px var(--c)}
.tbtn.on{color:var(--c);background:color-mix(in srgb,var(--c) 10%,var(--card))}
/* secciones */
section.blk{padding:84px 0 0}
.sh{display:flex;align-items:flex-end;justify-content:space-between;gap:20px;margin-bottom:28px;flex-wrap:wrap}
.sh .n{font:600 12px var(--mono);letter-spacing:.2em;color:var(--burn);text-transform:uppercase}
h2{font:700 clamp(26px,3.4vw,38px)/1.1 var(--disp);letter-spacing:-.02em;margin:8px 0 0}
.sh p{margin:0;color:var(--mut);max-width:520px;font-size:14px}
.card{background:var(--card);border:1px solid var(--line);border-radius:var(--r)}
/* narrativa */
.why{display:grid;grid-template-columns:1.6fr 1fr;gap:20px}
.story{padding:34px 36px;font-size:16.5px;line-height:1.75;color:var(--ink2)}
.story p{margin:0 0 16px}.story p:first-child{font:500 21px/1.5 var(--disp);color:var(--ink)}.story strong{color:#ffb38a;font-weight:600}
.usos{padding:28px}.usos h3{font:600 12px var(--mono);letter-spacing:.18em;text-transform:uppercase;color:var(--mut);margin:0 0 16px}
.uso{display:flex;gap:14px;align-items:flex-start;padding:14px 0;border-top:1px solid var(--line);font-size:14.5px;color:var(--ink)}
.uso:first-of-type{border-top:0}.uso span{flex:none;width:30px;height:30px;border-radius:9px;display:grid;place-items:center;background:rgba(193,68,14,.15);color:var(--burn);font:600 13px var(--mono)}
/* fortalezas y debilidades */
.sw{display:grid;grid-template-columns:repeat(3,1fr);gap:16px}
.swc{padding:26px}.swc h3{display:flex;align-items:center;gap:10px;font:600 17px var(--disp);margin:0 0 4px}
.swc .s{font-size:12.5px;color:var(--mut);margin-bottom:16px}
.swc ul{list-style:none;margin:0;padding:0}.swc li{padding:11px 0;border-top:1px solid var(--line);font-size:14px;display:flex;justify-content:space-between;gap:12px}
.swc li span{color:var(--ink2)}.swc li b{font:600 13px var(--mono);white-space:nowrap}
.swc.ok h3 i,.swc.mid h3 i,.swc.bad h3 i{width:10px;height:10px;border-radius:50%}
.swc.ok h3 i{background:var(--ok);box-shadow:0 0 10px var(--ok)}.swc.mid h3 i{background:var(--mid);box-shadow:0 0 10px var(--mid)}.swc.bad h3 i{background:var(--bad);box-shadow:0 0 10px var(--bad)}
.swc.ok li b{color:var(--ok)}.swc.mid li b{color:var(--mid)}.swc.bad li b{color:var(--bad)}
.swc li.more{display:none}.swc.all li.more{display:flex}.moreb{margin-top:14px;width:100%}
.empty{color:var(--mut2);font-size:13.5px;padding:10px 0}
/* puntuación por sección */
.scores{display:grid;grid-template-columns:repeat(2,1fr);gap:14px}
.sc{display:grid;grid-template-columns:1fr auto;align-items:center;gap:6px 18px;padding:22px 24px;cursor:pointer;transition:.25s;text-decoration:none}
.sc:hover{border-color:rgba(224,120,47,.45);background:var(--card2)}
.sc .t{font:600 16px var(--disp)}.sc .m{font-size:12.5px;color:var(--mut)}
.sc .num{grid-row:span 2;font:700 44px/1 var(--disp);letter-spacing:-.03em;text-align:right}
.sc .bar{grid-column:1/-1;height:8px;margin-top:8px}
/* filtros + acordeón */
.tools{display:flex;flex-wrap:wrap;gap:10px;align-items:center;margin-bottom:18px}
.tools .sp{flex:1}
.btn{font:500 13px var(--body);color:var(--ink2);background:var(--card);border:1px solid var(--line);padding:9px 14px;border-radius:10px;cursor:pointer;transition:.2s}
.btn:hover{color:var(--ink);border-color:var(--line2)}
.search{font:400 14px var(--body);color:var(--ink);background:var(--card);border:1px solid var(--line);padding:9px 14px;border-radius:10px;min-width:220px;outline:none}
.search:focus{border-color:var(--burn)}
.acc{margin-bottom:12px;overflow:hidden;transition:border-color .2s}
.acc.open{border-color:var(--line2)}
.acch{width:100%;display:grid;grid-template-columns:auto 1fr auto auto 26px;align-items:center;gap:18px;padding:20px 24px;background:none;border:0;color:inherit;cursor:pointer;text-align:left;font:inherit}
.acch:hover{background:rgba(255,255,255,.02)}
.acch .ix{font:600 12px var(--mono);color:var(--mut2)}
.acch .t{font:600 19px var(--disp)}.acch .t small{display:block;font:400 12.5px var(--body);color:var(--mut);margin-top:2px}
.minis{display:flex;gap:6px}.mini{font:600 11.5px var(--mono);padding:4px 8px;border-radius:7px}
.mini.ok{color:var(--ok);background:var(--okbg)}.mini.mid{color:var(--mid);background:var(--midbg)}.mini.bad{color:var(--bad);background:var(--badbg)}.mini.na{color:var(--na);background:var(--nabg)}
.acch .scw{width:150px}.acch .scw b{font:700 22px var(--disp);display:block;text-align:right}.acch .scw .bar{margin-top:6px}
.chev{width:26px;height:26px;border-radius:8px;border:1px solid var(--line2);display:grid;place-items:center;transition:transform .3s}
.chev:before{content:"";width:7px;height:7px;border-right:2px solid var(--ink2);border-bottom:2px solid var(--ink2);transform:translateY(-2px) rotate(45deg)}
.acc.open .chev{transform:rotate(180deg)}
.accb{display:grid;grid-template-rows:0fr;transition:grid-template-rows .4s ease}
.acc.open .accb{grid-template-rows:1fr}
.accb>div{overflow:hidden}
/* filas de variables */
.vh,.vr{display:grid;grid-template-columns:2.1fr 1.3fr 1.4fr 118px 56px 64px 44px;gap:16px;align-items:start;padding:15px 24px}
.vh{font:600 10.5px var(--mono);letter-spacing:.14em;text-transform:uppercase;color:var(--mut2);border-top:1px solid var(--line);padding-top:12px;padding-bottom:12px}
.vr{border-top:1px solid var(--line);transition:background .2s}
.vr:hover{background:rgba(255,255,255,.025)}
.vr .vn{font-weight:600;font-size:14.5px;display:flex;gap:8px;align-items:flex-start}
.vr .val{font:600 14px var(--mono);color:var(--ink)}.vr .rng{font:400 13px var(--mono);color:var(--ink2)}
.vr .note{font:400 12px/1.5 var(--body);color:var(--mut);margin-top:5px}
.vr .sco{font:600 14px var(--mono);text-align:right}.vr .rf{font:500 12px var(--mono);color:var(--mut);text-align:center;white-space:nowrap}
.vr .rf b{color:var(--ink2);font-weight:600}
.src{grid-column:1/-1;display:none;font-size:12.5px;color:var(--mut);background:var(--bg2);border:1px solid var(--line);border-radius:10px;padding:12px 14px;margin-top:-4px;line-height:1.55}
.src b{color:var(--ink2);font-weight:600}
.vr.showsrc .src{display:block}
.srcb{background:none;border:1px solid var(--line2);color:var(--mut);width:34px;height:28px;border-radius:8px;cursor:pointer;font:600 11px var(--mono);transition:.2s}
.srcb:hover,.vr.showsrc .srcb{color:var(--burn);border-color:var(--burn)}
.vr.hide,.acc.hide{display:none}
/* pills de estado */
.pill{display:inline-flex;align-items:center;gap:7px;font:600 12px var(--body);padding:5px 11px 5px 9px;border-radius:999px;white-space:nowrap}
.pill:before{content:"";width:7px;height:7px;border-radius:50%;background:currentColor}
.pill.dentro{color:var(--ok);background:var(--okbg)}.pill.cercano{color:var(--mid);background:var(--midbg)}
.pill.fuera{color:var(--bad);background:var(--badbg)}.pill.sin,.pill.cualitativo{color:var(--na);background:var(--nabg)}
/* tooltip */
.tip{flex:none;width:18px;height:18px;border-radius:50%;border:1px solid var(--line2);color:var(--mut);display:inline-grid;place-items:center;font:600 10.5px var(--mono);cursor:help;margin-top:2px;transition:.2s;background:none;padding:0}
.tip:hover,.tip:focus{color:var(--burn);border-color:var(--burn);outline:none}
#tt{position:fixed;z-index:60;max-width:340px;background:#1b1f27;border:1px solid var(--line2);color:var(--ink2);font:400 13px/1.55 var(--body);padding:12px 14px;border-radius:12px;
  box-shadow:0 20px 50px -10px #000;pointer-events:none;opacity:0;transform:translateY(4px);transition:opacity .15s,transform .15s}
#tt.on{opacity:1;transform:none}#tt b{display:block;color:var(--burn);font:600 10.5px var(--mono);letter-spacing:.14em;text-transform:uppercase;margin-bottom:5px}
/* límites, misión, pendientes */
.grid3{display:grid;grid-template-columns:repeat(3,1fr);gap:14px}
.lim{padding:22px;transition:.25s}.lim:hover{border-color:var(--line2);transform:translateY(-2px)}
.lim .t{font:600 16px var(--disp);display:flex;justify-content:space-between;gap:10px}
.lim .mv{font:700 24px/1.2 var(--disp);color:var(--burn);margin:12px 0 4px}
.lim .ml{font:600 10.5px var(--mono);letter-spacing:.14em;text-transform:uppercase;color:var(--mut2)}
.lim .sv{font:500 13px var(--mono);color:var(--ink2);margin-top:10px}
.lim p{font-size:12.5px;color:var(--mut);margin:12px 0 0;line-height:1.55}
.list{overflow:hidden}
.lr{display:grid;grid-template-columns:1.4fr 1fr 1.6fr;gap:18px;padding:14px 24px;border-top:1px solid var(--line);font-size:14px;transition:background .2s}
.lr:first-child{border-top:0}.lr:hover{background:rgba(255,255,255,.025)}
.lr .k{font-weight:600}.lr .m{font:400 13px var(--mono);color:var(--ink2)}.lr .s{font-size:12.5px;color:var(--mut)}
.lr.hd{font:600 10.5px var(--mono);letter-spacing:.14em;text-transform:uppercase;color:var(--mut2);background:var(--bg2)}
.sub2{display:flex;align-items:center;gap:14px;margin:48px 0 16px}
.sub2 h3{font:600 20px var(--disp);margin:0}.sub2 span{font:600 11px var(--mono);color:var(--mut);background:var(--card2);border:1px solid var(--line);padding:4px 9px;border-radius:999px}
/* spaceport */
.sp-head{display:grid;grid-template-columns:auto 1fr;gap:32px;align-items:center;padding:28px 30px;margin-bottom:14px}
.sp-head .num{font:700 64px/1 var(--disp);letter-spacing:-.04em}.sp-head .num small{display:block;font:600 11px var(--mono);letter-spacing:.18em;color:var(--mut);text-transform:uppercase;margin-top:8px}
.sp-head p{margin:0 0 12px;color:var(--ink2);font-size:14px}
.spr{display:grid;grid-template-columns:1.3fr .9fr 1.3fr 118px;gap:16px;padding:16px 24px;border-top:1px solid var(--line);align-items:start;transition:background .2s}
.spr:hover{background:rgba(255,255,255,.025)}.spr:first-child{border-top:0}
.spr .k{font-weight:600;font-size:14.5px;display:flex;gap:8px}.spr .m{font:600 13.5px var(--mono)}.spr .c{font:400 13px var(--mono);color:var(--ink2)}
.spr .alt{grid-column:1/-1;font-size:12.5px;color:var(--mut);display:none;padding:10px 14px;background:var(--bg2);border:1px solid var(--line);border-radius:10px}
.spr .alt b{color:var(--burn);font-weight:600}.spr.showsrc .alt{display:block}
/* atributos + notas */
.kv{display:grid;grid-template-columns:repeat(auto-fill,minmax(220px,1fr));gap:12px}
.kv div{padding:16px 18px}.kv span{display:block;font:600 10.5px var(--mono);letter-spacing:.14em;text-transform:uppercase;color:var(--mut2)}.kv b{display:block;font:500 14px var(--mono);margin-top:6px;word-break:break-word}
.notes{padding:8px 30px}.notes li{padding:13px 0;border-top:1px solid var(--line);color:var(--ink2);font-size:14px}.notes li:first-child{border-top:0}
.notes ul{list-style:none;margin:0;padding:0}
/* footer */
footer{margin-top:110px;border-top:1px solid var(--line);padding:38px 0 56px}
footer .wrap{display:flex;justify-content:space-between;gap:24px;flex-wrap:wrap;align-items:center}
footer .b{font:600 13px var(--mono);letter-spacing:.16em;text-transform:uppercase;color:var(--ink2);display:flex;align-items:center;gap:10px}
footer .b i{width:8px;height:8px;border-radius:50%;background:var(--mars)}
footer p{margin:0;font-size:12.5px;color:var(--mut)}
/* volver arriba */
#up{position:fixed;right:24px;bottom:24px;z-index:50;width:50px;height:50px;border-radius:50%;border:1px solid rgba(224,120,47,.5);background:rgba(18,21,27,.85);backdrop-filter:blur(10px);
  color:var(--burn);cursor:pointer;display:grid;place-items:center;opacity:0;transform:translateY(16px);pointer-events:none;transition:.3s}
#up.on{opacity:1;transform:none;pointer-events:auto}#up:hover{background:var(--mars);color:#fff;box-shadow:0 10px 30px -6px var(--glow)}
#up:before{content:"";width:10px;height:10px;border-left:2px solid currentColor;border-top:2px solid currentColor;transform:translateY(3px) rotate(45deg)}
/* aparición al hacer scroll */
.rv{opacity:0;transform:translateY(22px);transition:opacity .8s ease,transform .8s ease}.rv.in{opacity:1;transform:none}
@media (prefers-reduced-motion:reduce){*{transition:none!important;animation:none!important}.rv{opacity:1;transform:none}}
/* responsive */
@media (max-width:980px){
  .hgrid,.why{grid-template-columns:1fr}.gauge{justify-self:center;width:260px}
  .kpis{grid-template-columns:repeat(2,1fr)}.sw,.grid3{grid-template-columns:1fr}.scores{grid-template-columns:1fr}
  .vh{display:none}
  .vr{grid-template-columns:1fr auto;gap:8px 14px}
  .vr .vn{grid-column:1/-1}.vr .c-val,.vr .c-rng{grid-column:1/-1}
  .vr .c-val:before,.vr .c-rng:before{display:block;font:600 10px var(--mono);letter-spacing:.14em;text-transform:uppercase;color:var(--mut2);margin-bottom:2px}
  .vr .c-val:before{content:"Valor del sitio"}.vr .c-rng:before{content:"Ventana Marte"}
  .vr .sco,.vr .rf{text-align:left}
  .acch{grid-template-columns:1fr auto 26px;gap:12px}.acch .ix,.acch .minis{display:none}.acch .scw{width:90px}
  .spr,.lr{grid-template-columns:1fr}.sp-head{grid-template-columns:1fr}
}
@media (max-width:560px){.wrap{padding:0 18px}.hero{padding-top:56px}.kpis{gap:10px}.kpi{padding:16px}.kpi .v{font-size:28px}.story{padding:24px}.nav{display:none}}
@media print{.top,#up,.tools,.tally{display:none}body{background:#fff;color:#000}.accb{grid-template-rows:1fr!important}.rv{opacity:1;transform:none}}
"""

FICHA_JS = r"""
(function(){
  const $=(s,c=document)=>c.querySelector(s), $$=(s,c=document)=>[...c.querySelectorAll(s)];
  // barras y aparición al entrar en pantalla
  const io=new IntersectionObserver(es=>es.forEach(e=>{if(!e.isIntersecting)return;const t=e.target;
    if(t.classList.contains('rv'))t.classList.add('in');
    $$('.bar>i[data-w]',t).forEach(b=>b.style.width=b.dataset.w+'%');
    if(t.matches('.bar>i[data-w]'))t.style.width=t.dataset.w+'%';
    io.unobserve(t);}),{threshold:.12});
  $$('.rv, .bar>i[data-w]').forEach(el=>io.observe(el));
  // indicador circular + contador
  const g=$('.gauge .val');if(g){const C=2*Math.PI*54,v=+g.dataset.v||0;g.style.strokeDasharray=C;g.style.strokeDashoffset=C;
    requestAnimationFrame(()=>setTimeout(()=>g.style.strokeDashoffset=C*(1-v),150));}
  $$('[data-count]').forEach(el=>{const end=parseFloat(el.dataset.count),dec=+(el.dataset.dec||0),t0=performance.now(),D=1600;
    const step=t=>{const p=Math.min(1,(t-t0)/D),k=1-Math.pow(1-p,3);el.textContent=(end*k).toFixed(dec);if(p<1)requestAnimationFrame(step)};requestAnimationFrame(step);});
  // acordeones
  const setOpen=(a,o)=>{a.classList.toggle('open',o);$('.acch',a).setAttribute('aria-expanded',o)};
  $$('.acc').forEach(a=>$('.acch',a).addEventListener('click',()=>{const o=!a.classList.contains('open');setOpen(a,o);if(o)$$('.bar>i[data-w]',a).forEach(b=>b.style.width=b.dataset.w+'%')}));
  $('#openall')&&$('#openall').addEventListener('click',()=>$$('.acc:not(.hide)').forEach(a=>setOpen(a,true)));
  $('#closeall')&&$('#closeall').addEventListener('click',()=>$$('.acc').forEach(a=>setOpen(a,false)));
  $$('.sc[data-go]').forEach(s=>s.addEventListener('click',e=>{e.preventDefault();const a=document.getElementById(s.dataset.go);if(a){setOpen(a,true);a.scrollIntoView({behavior:'smooth'})}}));
  $$('.moreb').forEach(b=>b.addEventListener('click',()=>{const c=b.closest('.swc');c.classList.toggle('all');b.textContent=c.classList.contains('all')?'Ver menos':'Ver '+b.dataset.n+' más'}));
  // fuentes por fila
  $$('.srcb').forEach(b=>b.addEventListener('click',()=>{const r=b.closest('.vr,.spr');r.classList.toggle('showsrc');b.setAttribute('aria-expanded',r.classList.contains('showsrc'))}));
  // filtros por estado + búsqueda
  let filtro=null;const q=$('#q');
  function aplicar(){const txt=(q&&q.value||'').trim().toLowerCase();
    $$('.acc').forEach(a=>{let n=0;$$('.vr',a).forEach(r=>{const ok=(!filtro||r.dataset.estado===filtro)&&(!txt||r.textContent.toLowerCase().includes(txt));r.classList.toggle('hide',!ok);if(ok)n++});
      a.classList.toggle('hide',n===0);if((filtro||txt)&&n)setOpen(a,true);});
    $$('.tbtn').forEach(b=>b.classList.toggle('on',b.dataset.f===filtro));
    const lab=$('#flab');if(lab)lab.textContent=filtro?('Filtro: '+filtro):'';}
  $$('.tbtn').forEach(b=>b.addEventListener('click',()=>{filtro=filtro===b.dataset.f?null:b.dataset.f;aplicar();if(b.closest('.hero'))document.getElementById('variables').scrollIntoView({behavior:'smooth'})}));
  q&&q.addEventListener('input',aplicar);
  $('#clearf')&&$('#clearf').addEventListener('click',()=>{filtro=null;if(q)q.value='';aplicar();$$('.acc').forEach(a=>setOpen(a,false))});
  // tooltips
  const tt=$('#tt');let cur=null;
  function show(el){cur=el;tt.innerHTML='<b>'+(el.dataset.h||'Por qué importa')+'</b>'+el.dataset.tip.replace(/</g,'&lt;');tt.classList.add('on');
    const r=el.getBoundingClientRect(),w=tt.offsetWidth,h=tt.offsetHeight;let x=r.left+r.width/2-w/2,y=r.top-h-10;
    x=Math.max(10,Math.min(x,innerWidth-w-10));if(y<64)y=r.bottom+10;tt.style.left=x+'px';tt.style.top=y+'px';}
  function hide(){cur=null;tt.classList.remove('on')}
  $$('.tip').forEach(el=>{el.addEventListener('mouseenter',()=>show(el));el.addEventListener('mouseleave',hide);el.addEventListener('focus',()=>show(el));el.addEventListener('blur',hide);
    el.addEventListener('click',e=>{e.stopPropagation();cur===el?hide():show(el)})});
  addEventListener('scroll',()=>{if(cur)hide()},{passive:true});document.addEventListener('click',hide);
  // volver arriba + navegación activa
  const up=$('#up');up.addEventListener('click',()=>scrollTo({top:0,behavior:'smooth'}));
  const secs=$$('section[id]'),links=$$('.nav a');
  addEventListener('scroll',()=>{up.classList.toggle('on',scrollY>700);let id='';secs.forEach(s=>{if(s.getBoundingClientRect().top<140)id=s.id});
    links.forEach(l=>l.classList.toggle('on',l.getAttribute('href')==='#'+id));},{passive:true});
})();
"""

ESTADO_CLS = {"dentro": "dentro", "cercano": "cercano", "fuera": "fuera", "sin dato": "sin", "cualitativo": "cualitativo"}
ESTADO_TXT = {"dentro": "Dentro", "cercano": "Cercano", "fuera": "Fuera", "sin dato": "Sin dato", "cualitativo": "Cualitativo"}


def _pill(estado):
    return f"<span class='pill {ESTADO_CLS.get(estado, 'sin')}'>{ESTADO_TXT.get(estado, estado)}</span>"


def _barcls(s):
    if s is None:
        return ""
    return "ok" if s >= 0.8 else "mid" if s >= 0.5 else "bad"


def _bar(s, cls=None):
    w = 0 if s is None else max(0, min(1, s)) * 100
    return f"<div class='bar {cls if cls is not None else _barcls(s)}'><i data-w='{w:.0f}'></i></div>"


def _tip(texto, titulo="Por qué importa"):
    if not texto:
        return ""
    return f"<button type='button' class='tip' aria-label='{e(titulo)}' data-h='{e(titulo)}' data-tip='{e(texto)}'>i</button>"


def ficha_html(f, notas):
    g = f["geo"] or {}
    i = f["indice"]
    vars_eval = [r for r in f["variables"] if r["tipo"] in ("numerica", "categorica")]
    cuenta = {k: sum(1 for r in vars_eval if r["estado"] == k) for k in ("dentro", "cercano", "fuera", "sin dato")}
    limites = [r for r in f["variables"] if r["tipo"] == "limite"]
    falt = [r for r in vars_eval if r["score"] is None]
    mision = [r for r in f["variables"] if r["tipo"] in ("mision", "norma")]
    fuertes = sorted([r for r in vars_eval if r["estado"] == "dentro"], key=lambda r: -peso(r))
    cerca = sorted([r for r in vars_eval if r["estado"] == "cercano"], key=lambda r: -peso(r))
    debiles = sorted([r for r in vars_eval if r["estado"] == "fuera"], key=lambda r: -peso(r))
    secciones = []
    for r in vars_eval:
        if r["seccion"] not in secciones:
            secciones.append(r["seccion"])
    sid = {s: "sec-" + slug(s) for s in secciones}

    H = ["<!DOCTYPE html><html lang='es'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>",
         f"<title>{e(f['nombre'])} · Ficha científica de análogo marciano</title>",
         f"<meta name='description' content='Ficha científica del análogo Tierra–Marte {e(f['nombre'])}: índice de analogía {fmt(i)}.'>",
         f"<style>{FICHA_CSS}</style></head><body>"]
    # barra superior
    nav = [("resumen", "Resumen"), ("balance", "Balance"), ("secciones", "Secciones"), ("variables", "Variables"),
           ("limites", "Límites")]
    if f["spaceport"]:
        nav.append(("spaceport", "Spaceport"))
    nav.append(("notas", "Notas"))
    H.append(f"<header class='top'><div class='wrap'><div class='brand'><i></i>Análogo Marte · {e(f['id'])}</div><nav class='nav'>"
             + "".join(f"<a href='#{a}'>{t}</a>" for a, t in nav) + "</nav></div></header>")

    # hero
    H.append("<section class='hero' id='top'><div class='wrap'><div class='hgrid'><div>"
             "<div class='kick'>Ficha científica · Análogo Tierra–Marte</div>"
             f"<h1>{e(f['nombre'])}</h1><div class='chips'><span class='chip'>ID {e(f['id'])}</span>")
    if g:
        H.append(f"<span class='chip'>{g['lat']:.5f}°, {g['lon']:.5f}°</span><span class='chip'>Extensión {fmt(g['extension_km'])} km</span>"
                 f"<span class='chip'>Área {fmt(g['area_km2'])} km²</span>")
    H.append(f"</div><div class='cls'><i></i>{e(f['clase'])}</div></div>")
    gv = 0 if i is None else max(0, min(1, i))
    H.append("<div class='gauge'><svg viewBox='0 0 120 120'><defs><linearGradient id='gg' x1='0' y1='0' x2='1' y2='1'>"
             "<stop offset='0' stop-color='#e0782f'/><stop offset='1' stop-color='#c1440e'/></linearGradient></defs>"
             "<circle class='trk' cx='60' cy='60' r='54' fill='none' stroke-width='7'/>"
             f"<circle class='val' data-v='{gv:.4f}' cx='60' cy='60' r='54' fill='none' stroke-width='7'/></svg>"
             f"<div class='ctr'><div class='big' data-count='{gv:.3f}' data-dec='2'>{fmt(i)}</div>"
             "<div class='lab'>Índice de analogía</div><div class='sub'>escala 0 – 1</div></div></div></div>")
    # KPIs
    k = [("Índice de analogía", f"<span data-count='{gv:.3f}' data-dec='3'>{fmt(i)}</span>", "Ponderado por relevancia × fidelidad", i),
         ("Cobertura de datos", f"<span data-count='{f['cobertura']*100:.0f}'>{f['cobertura']*100:.0f}</span><small>%</small>",
          f"{f['n_evaluadas']} de {f['n_evaluables']} variables evaluables", f["cobertura"])]
    if f["indice_spaceport"] is not None:
        k.append(("Aptitud spaceport", f"<span data-count='{f['indice_spaceport']:.3f}' data-dec='3'>{fmt(f['indice_spaceport'])}</span>",
                  "Criterios terrestres del equipo", f["indice_spaceport"]))
    k.append(("Variables dentro", f"<span data-count='{cuenta['dentro']}'>{cuenta['dentro']}</span><small>/ {len(vars_eval)}</small>",
              "Dentro de la ventana marciana", cuenta["dentro"] / len(vars_eval) if vars_eval else 0))
    if g and len(k) < 4:
        k.append(("Área", f"{fmt(g['area_km2'])}<small>km²</small>", f"Extensión {fmt(g['extension_km'])} km", None))
    H.append("<div class='kpis'>")
    for t, v, d, b in k:
        H.append(f"<div class='kpi'><div class='l'>{t}</div><div class='v'>{v}</div><div class='d'>{e(d)}</div>"
                 + (_bar(b, "") if b is not None else "") + "</div>")
    H.append("</div><div class='tally'>")
    for est in ("dentro", "cercano", "fuera", "sin dato"):
        H.append(f"<button type='button' class='tbtn' data-f='{est}'><span class='dot'></span><b>{cuenta[est]}</b>{ESTADO_TXT[est]}</button>")
    H.append("</div></div></section>")

    # resumen
    H.append("<main class='wrap'><section class='blk rv' id='resumen'><div class='sh'><div><div class='n'>01 · Resumen</div>"
             "<h2>¿Por qué es análogo de Marte?</h2></div></div><div class='why'><div class='card story'>")
    for p in f["narrativa"]:
        H.append("<p>" + re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", e(p)) + "</p>")
    H.append("</div><div class='card usos'><h3>Usos recomendados del análogo</h3>")
    if f["usos_recomendados"]:
        for n_, u in enumerate(f["usos_recomendados"], 1):
            H.append(f"<div class='uso'><span>{n_:02d}</span><div>{e(u)}</div></div>")
    else:
        H.append("<div class='empty'>Ninguna sección alcanza aún el umbral (≥ 0.7 con al menos 2 variables).</div>")
    H.append("</div></div></section>")

    # balance
    def lista(rs, vacio):
        if not rs:
            return f"<div class='empty'>{vacio}</div>"
        li = "".join(f"<li{' class=\'more\'' if j >= 6 else ''}><span>{e(r['variable'])}</span><b>{e(valor_txt(r))}</b></li>"
                     for j, r in enumerate(rs))
        extra = f"<button type='button' class='btn moreb' data-n='{len(rs) - 6}'>Ver {len(rs) - 6} más</button>" if len(rs) > 6 else ""
        return "<ul>" + li + "</ul>" + extra
    H.append("<section class='blk rv' id='balance'><div class='sh'><div><div class='n'>02 · Balance</div><h2>Fortalezas y debilidades</h2></div>"
             "<p>Variables con dato ordenadas por peso (relevancia × fidelidad).</p></div><div class='sw'>"
             f"<div class='card swc ok'><h3><i></i>Fortalezas</h3><div class='s'>{len(fuertes)} variables dentro de la ventana marciana</div>{lista(fuertes, 'Sin variables dentro de la ventana.')}</div>"
             f"<div class='card swc mid'><h3><i></i>Se acerca</h3><div class='s'>{len(cerca)} variables cercanas, sin igualar</div>{lista(cerca, 'Sin variables cercanas.')}</div>"
             f"<div class='card swc bad'><h3><i></i>Debilidades</h3><div class='s'>{len(debiles)} variables fuera de la ventana</div>{lista(debiles, 'Sin variables fuera de la ventana.')}</div>"
             "</div></section>")

    # puntuación por sección
    H.append("<section class='blk rv' id='secciones'><div class='sh'><div><div class='n'>03 · Desempeño</div><h2>Puntuación por sección</h2></div>"
             "<p>Media ponderada por relevancia × fidelidad de las variables con dato. Haz clic en una sección para ver sus variables.</p></div><div class='scores'>")
    for s, (sc, n) in f["secciones"].items():
        H.append(f"<a class='card sc' href='#{sid.get(s, '')}' data-go='{sid.get(s, '')}'><div class='t'>{e(s)}</div>"
                 f"<div class='num'>{fmt(sc)}</div><div class='m'>{n} variable{'s' if n != 1 else ''} con dato</div>{_bar(sc)}</a>")
    H.append("</div></section>")

    # variables
    H.append("<section class='blk' id='variables'><div class='sh'><div><div class='n'>04 · Propiedades</div>"
             "<h2>Propiedades del sitio frente a la ventana marciana</h2></div>"
             "<p>R = relevancia para el diseño (1–5) · F = fidelidad del análogo terrestre (1–5). Pasa el cursor sobre «i» para ver por qué importa cada variable.</p></div>"
             "<div class='tools'>")
    for est in ("dentro", "cercano", "fuera", "sin dato"):
        H.append(f"<button type='button' class='tbtn' data-f='{est}'><span class='dot'></span><b>{cuenta[est]}</b>{ESTADO_TXT[est]}</button>")
    H.append("<span class='sp'></span><input id='q' class='search' type='search' placeholder='Buscar variable…' aria-label='Buscar variable'>"
             "<button type='button' class='btn' id='openall'>Expandir todo</button><button type='button' class='btn' id='closeall'>Contraer</button>"
             "<button type='button' class='btn' id='clearf'>Limpiar</button></div>")
    for n_, s in enumerate(secciones, 1):
        rs = [r for r in vars_eval if r["seccion"] == s]
        sc = f["secciones"].get(s, (None, 0))[0]
        c = {k2: sum(1 for r in rs if r["estado"] == k2) for k2 in ("dentro", "cercano", "fuera", "sin dato")}
        minis = "".join(f"<span class='mini {cl}'>{c[k2]} {lb}</span>" for k2, cl, lb in
                        (("dentro", "ok", "✓"), ("cercano", "mid", "≈"), ("fuera", "bad", "✕"), ("sin dato", "na", "—")) if c[k2])
        H.append(f"<div class='card acc' id='{sid[s]}'><button type='button' class='acch' aria-expanded='false'>"
                 f"<span class='ix'>{n_:02d}</span><span class='t'>{e(s)}<small>{len(rs)} variables · {c['dentro']} dentro de la ventana marciana</small></span>"
                 f"<span class='minis'>{minis}</span><span class='scw'><b>{fmt(sc)}</b>{_bar(sc)}</span><span class='chev'></span></button>"
                 "<div class='accb'><div><div class='vh'><span>Variable</span><span>Valor del sitio</span><span>Ventana Marte</span><span>Estado</span>"
                 "<span style='text-align:right'>Punt.</span><span style='text-align:center'>R · F</span><span></span></div>")
        for r in rs:
            nota = " ".join(x for x in (r["detalle"], r["observacion"]) if x)
            H.append(f"<div class='vr' data-estado='{e(r['estado'])}'>"
                     f"<div class='vn'><span>{e(r['variable'])}</span>{_tip(r['por_que'])}</div>"
                     f"<div class='c-val'><div class='val'>{e(valor_txt(r))}</div>" + (f"<div class='note'>{e(nota)}</div>" if nota else "") + "</div>"
                     f"<div class='c-rng rng'>{e(rango_txt(r))}</div><div>{_pill(r['estado'])}</div>"
                     f"<div class='sco'>{fmt(r['score'])}</div><div class='rf'><b>{fmt(r['relevancia'])}</b> · <b>{fmt(r['fidelidad'])}</b></div>"
                     "<div><button type='button' class='srcb' aria-expanded='false' title='Ver fuentes NASA'>SRC</button></div>"
                     f"<div class='src'><b>Fuente Marte (NASA):</b> {e(r['fuente_marte'])}<br><b>Fuente análogo:</b> {e(r['fuente_analogo'])}"
                     + (f"<br><b>Columna de la capa:</b> <code>{e(r['columna'])}</code>" if r["columna"] else "") + "</div></div>")
        H.append("</div></div></div>")
    H.append("</section>")

    # límites
    H.append("<section class='blk rv' id='limites'><div class='sh'><div><div class='n'>05 · Límites</div><h2>Lo que la Tierra no puede imitar</h2></div>"
             "<p>Variables relevantes para el diseño en Marte, pero con fidelidad 1. Se declaran como límites de la analogía y se resuelven con cámaras, simulación o laboratorio.</p></div><div class='grid3'>")
    for r in limites:
        H.append(f"<div class='card lim'><div class='t'><span>{e(r['variable'])}</span>{_tip(r['por_que'])}</div>"
                 f"<div class='mv'>{e(rango_txt(r))}</div><div class='ml'>Valor en Marte · R {fmt(r['relevancia'])}</div>"
                 f"<div class='sv'>Sitio: {e(valor_txt(r))}</div><p>{e(r['observacion'] or r['fuente_analogo'])}</p>"
                 f"<p><b style='color:var(--ink2)'>Fuente:</b> {e(r['fuente_marte'])}</p></div>")
    H.append("</div>")
    if falt:
        H.append(f"<div class='sub2'><h3>Pendiente de medir</h3><span>{len(falt)} variables</span></div><div class='card list'>"
                 "<div class='lr hd'><span>Variable</span><span>Ventana Marte</span><span>Cómo obtenerla</span></div>")
        for r in falt:
            H.append(f"<div class='lr'><span class='k'>{e(r['variable'])}</span><span class='m'>{e(rango_txt(r))}</span><span class='s'>{e(r['fuente_analogo'])}</span></div>")
        H.append("</div>")
    H.append(f"<div class='sub2'><h3>Parámetros de misión</h3><span>no dependen del sitio</span></div><div class='card list'>"
             "<div class='lr hd'><span>Parámetro</span><span>Valor de referencia</span><span>Fuente</span></div>")
    for r in mision:
        H.append(f"<div class='lr'><span class='k'>{e(r['variable'])}</span><span class='m'>{e(rango_txt(r))}</span><span class='s'>{e(r['fuente_marte'])}"
                 + (f" · {e(r['observacion'] or r['fuente_analogo'])}" if (r['observacion'] or r['fuente_analogo']) else "") + "</span></div>")
    H.append("</div></section>")

    # spaceport
    if f["spaceport"]:
        sp = f["spaceport"]
        sc_c = {k2: sum(1 for r in sp if r["estado"] == k2) for k2 in ("dentro", "cercano", "fuera", "sin dato", "cualitativo")}
        H.append("<section class='blk rv' id='spaceport'><div class='sh'><div><div class='n'>06 · Operaciones en Tierra</div><h2>Aptitud como spaceport terrestre</h2></div>"
                 "<p>Criterios del equipo (documento «Variables por definir»). Pulsa ALT para ver la alternativa de ingeniería.</p></div>"
                 f"<div class='card sp-head'><div class='num'>{fmt(f['indice_spaceport'])}<small>Índice spaceport</small></div><div>"
                 f"<p>{sc_c['dentro']} criterios cumplidos · {sc_c['cercano']} cercanos · {sc_c['fuera']} fuera · {sc_c['sin dato']} sin dato · {sc_c['cualitativo']} cualitativos</p>"
                 f"{_bar(f['indice_spaceport'])}</div></div><div class='card list'>"
                 "<div class='spr lr hd' style='display:grid'><span>Parámetro</span><span>Valor</span><span>Criterio</span><span>Estado</span></div>")
        for r in sp:
            crit = r["ideal"] or f"{fmt(r['min'])} a {fmt(r['max'])} {r['unidad']}"
            val = r["valor"] if isinstance(r["valor"], str) else (f"{fmt(r['valor'])} {r['unidad']}" if r["valor"] is not None else "—")
            alt = f"<b>Alternativa:</b> {e(r['alternativa'])}" if r["alternativa"] else ""
            H.append(f"<div class='spr'><div class='k'><span>{e(r['parametro'])}</span>{_tip(r['utilidad'], 'Utilidad')}</div>"
                     f"<div class='m'>{e(val)}</div><div class='c'>{e(crit)}</div><div style='display:flex;gap:8px;align-items:center'>{_pill(r['estado'])}"
                     + (f"<button type='button' class='srcb' aria-expanded='false' title='Alternativa de ingeniería'>ALT</button>" if alt else "")
                     + "</div>" + (f"<div class='alt'>{alt}</div>" if alt else "") + "</div>")
        H.append("</div></section>")

    # atributos extra + notas
    if f["atributos_extra"]:
        H.append("<section class='blk rv' id='atributos'><div class='sh'><div><div class='n'>Datos</div><h2>Otros atributos de la capa</h2></div></div><div class='kv'>"
                 + "".join(f"<div class='card'><span>{e(k2)}</span><b>{e(v)}</b></div>" for k2, v in f["atributos_extra"].items()) + "</div></section>")
    H.append("<section class='blk rv' id='notas'><div class='sh'><div><div class='n'>Método</div><h2>Notas metodológicas</h2></div></div><div class='card notes'><ul>"
             + "".join(f"<li>{e(n_)}</li>" for n_ in notas)
             + "<li>Puntuación por variable: 1 dentro de la ventana marciana; fuera, decae linealmente con la distancia (tolerancia = ancho de la ventana o la definida en mapeo_columnas.csv). Categóricas: 1 si la clase observada coincide con una clase marciana.</li>"
             "<li>Índice de analogía = 30 % precipitación + 30 % NDVI + 40 % media de las demás variables ponderada por relevancia × fidelidad. Las variables de fidelidad 1 no puntúan.</li>"
             "</ul></div></section></main>")
    H.append("<footer><div class='wrap'><div class='b'><i></i>Análogos México · Dynamic Analog Engine</div>"
             f"<p>Ficha {e(f['id'])} · Generada por ficha_cientifica.py · Fuentes oficiales NASA · NASA Space Apps Challenge</p></div></footer>"
             "<button id='up' type='button' aria-label='Volver arriba'></button><div id='tt' role='tooltip'></div>"
             f"<script>{FICHA_JS}</script></body></html>")
    return "".join(H)


def kml(fichas):
    def coords(r):
        return " ".join(f"{x},{y},0" for x, y in r)
    P = ['<?xml version="1.0" encoding="UTF-8"?><kml xmlns="http://www.opengis.net/kml/2.2"><Document><name>Análogos Tierra–Marte</name>',
         '<Style id="a"><LineStyle><color>ff2b46b5</color><width>2</width></LineStyle><PolyStyle><color>662b46b5</color></PolyStyle></Style>']
    for f in fichas:
        top = [r for r in f["variables"] if r["score"] == 1 and r["tipo"] in ("numerica", "categorica")][:6]
        desc = (f"<b>Índice de analogía:</b> {fmt(f['indice'])} ({e(f['clase'])})<br><b>Cobertura:</b> {f['cobertura']*100:.0f} %<br>"
                + "<br>".join(f"✔ {e(r['variable'])}: {e(valor_txt(r))}" for r in top)
                + f"<br><a href='fichas/{e(f['id'])}.html'>Ver ficha completa</a>")
        P.append(f"<Placemark><name>{e(f['nombre'])}</name><styleUrl>#a</styleUrl><description><![CDATA[{desc}]]></description>")
        polys = [p for p in f["polys"] if p and len(p[0]) >= 3]
        if polys:
            P.append("<MultiGeometry>")
            for poly in polys:
                P.append("<Polygon><outerBoundaryIs><LinearRing><coordinates>" + coords(poly[0]) + "</coordinates></LinearRing></outerBoundaryIs>")
                for h in poly[1:]:
                    P.append("<innerBoundaryIs><LinearRing><coordinates>" + coords(h) + "</coordinates></LinearRing></innerBoundaryIs>")
                P.append("</Polygon>")
            P.append("</MultiGeometry>")
        elif f["geo"]:
            P.append(f"<Point><coordinates>{f['geo']['lon']},{f['geo']['lat']},0</coordinates></Point>")
        P.append("</Placemark>")
    P.append("</Document></kml>")
    return "".join(P)


def index_html(fichas):
    H = [f"<!DOCTYPE html><html lang='es'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>Fichas científicas — Análogos Tierra–Marte</title><style>{CSS}</style></head><body><div class='wrap'>",
         "<header><div class='kick'>NASA Space Apps · Análogos México</div><h1>Fichas científicas de sitios análogos a Marte</h1>",
         f"<div class='meta'><span>{len(fichas)} zonas evaluadas</span></div></header>",
         "<table><tr><th>#</th><th>Zona</th><th>Índice de analogía</th><th style='width:22%'></th><th>Clase</th><th>Cobertura</th><th>Spaceport</th></tr>"]
    for k, f in enumerate(fichas, 1):
        H.append(f"<tr><td>{k}</td><td><a href='fichas/{e(f['id'])}.html'><strong>{e(f['nombre'])}</strong></a></td><td>{fmt(f['indice'])}</td>"
                 f"<td><div class='bar'><i style='width:{(f['indice'] or 0)*100:.0f}%'></i></div></td><td>{e(f['clase'])}</td>"
                 f"<td>{f['cobertura']*100:.0f} %</td><td>{fmt(f['indice_spaceport'])}</td></tr>")
    H.append("</table></div></body></html>")
    return "".join(H)


def escribir_salidas(fichas, notas, out: Path):
    (out / "fichas").mkdir(parents=True, exist_ok=True)
    for f in fichas:
        (out / "fichas" / f"{f['id']}.html").write_text(ficha_html(f, notas), encoding="utf-8")
        (out / "fichas" / f"{f['id']}.md").write_text(ficha_md(f, notas), encoding="utf-8")
    limpio = [{k: v for k, v in f.items() if k != "polys"} for f in fichas]
    (out / "fichas.json").write_text(json.dumps({"notas": notas, "zonas": limpio}, ensure_ascii=False, indent=1, default=str), encoding="utf-8")
    with (out / "ranking.csv").open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        secs = sorted({s for f in fichas for s in f["secciones"]})
        w.writerow(["rank", "id", "nombre", "indice_analogia", "clase", "cobertura", "n_variables", "indice_spaceport", "lat", "lon", "area_km2"] + [f"sec_{s}" for s in secs] + ["usos_recomendados"])
        for k, f in enumerate(fichas, 1):
            g = f["geo"] or {}
            w.writerow([k, f["id"], f["nombre"], round(f["indice"], 4) if f["indice"] is not None else "", f["clase"],
                        round(f["cobertura"], 3), f["n_evaluadas"], round(f["indice_spaceport"], 3) if f["indice_spaceport"] is not None else "",
                        round(g.get("lat", 0), 6) if g else "", round(g.get("lon", 0), 6) if g else "", round(g.get("area_km2", 0), 3) if g else ""]
                       + [round(f["secciones"][s][0], 3) if s in f["secciones"] else "" for s in secs] + [" | ".join(f["usos_recomendados"])])
    feats = []
    for f in fichas:
        polys = [[[list(p) for p in r] for r in poly] for poly in f["polys"] if poly and len(poly[0]) >= 3]
        geom = {"type": "MultiPolygon", "coordinates": polys} if polys else (
            {"type": "Point", "coordinates": [f["geo"]["lon"], f["geo"]["lat"]]} if f["geo"] else None)
        feats.append({"type": "Feature", "geometry": geom, "properties": {
            "id": f["id"], "nombre": f["nombre"], "indice_analogia": f["indice"], "clase": f["clase"],
            "cobertura": f["cobertura"], "indice_spaceport": f["indice_spaceport"], "ficha": f"fichas/{f['id']}.html"}})
    (out / "zonas_resultado.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": feats}, ensure_ascii=False), encoding="utf-8")
    (out / "zonas_resultado.kml").write_text(kml(fichas), encoding="utf-8")
    (out / "index.html").write_text(index_html(fichas), encoding="utf-8")


# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description="Genera fichas científicas de sitios análogos Tierra–Marte.")
    ap.add_argument("--catalogo", required=True, type=Path, help="CSV de variables (Variables_Analogos_Tierra-Marte_corregido)")
    ap.add_argument("--zonas", type=Path, help="Capa poligonizada de QGIS (.geojson, .gpkg, .shp con geopandas, o .csv)")
    ap.add_argument("--capa", help="Nombre de la capa dentro del GeoPackage")
    ap.add_argument("--mapeo", type=Path, default=Path("mapeo_columnas.csv"))
    ap.add_argument("--spaceport", type=Path, default=Path("criterios_spaceport.csv"))
    ap.add_argument("--id-col", default="id")
    ap.add_argument("--nombre-col", default="nombre")
    ap.add_argument("--salida", type=Path, default=Path("salida"))
    ap.add_argument("--top", type=int, default=0, help="Solo las N zonas con mayor índice (0 = todas)")
    ap.add_argument("--crear-mapeo", action="store_true", help="Genera la plantilla mapeo_columnas.csv y termina")
    a = ap.parse_args(argv)

    variables, notas = leer_catalogo(a.catalogo)
    print(f"✓ Catálogo: {len(variables)} variables "
          f"({sum(v.tipo == 'numerica' for v in variables)} numéricas, {sum(v.tipo == 'categorica' for v in variables)} categóricas, "
          f"{sum(v.tipo == 'limite' for v in variables)} límites, {sum(v.tipo in ('mision', 'norma') for v in variables)} de misión)")
    if a.crear_mapeo:
        crear_mapeo(variables, a.mapeo)
        return
    if not a.zonas:
        ap.error("Falta --zonas")
    if a.mapeo.exists():
        aplicar_mapeo(variables, a.mapeo)
        print(f"✓ Mapeo aplicado: {a.mapeo}")
    else:
        print("⚠ No hay mapeo_columnas.csv: se usan los nombres de columna sugeridos.")
    criterios = leer_spaceport(a.spaceport)
    zonas = leer_zonas(a.zonas, a.capa)
    print(f"✓ Zonas leídas: {len(zonas)}")
    if zonas:
        cols = set(zonas[0]["props"])
        faltan = [v.columna for v in variables if v.tipo in ("numerica", "categorica") and v.columna and not buscar_col(zonas[0]["props"], v.columna) and v.columna not in cols
                  and not v.id.startswith(("semieje_mayor", "amplitud_termica"))]
        if faltan:
            print(f"  ℹ {len(faltan)} variables sin columna en la capa (irán a 'Pendiente de medir'): {', '.join(faltan[:12])}{'…' if len(faltan) > 12 else ''}")
    fichas = [construir_ficha(z, variables, criterios, a.id_col, a.nombre_col) for z in zonas]
    fichas.sort(key=lambda f: -(f["indice"] if f["indice"] is not None else -1))
    if a.top:
        fichas = fichas[: a.top]
    escribir_salidas(fichas, notas, a.salida)
    print(f"✓ {len(fichas)} fichas generadas en {a.salida}/  (abre {a.salida}/index.html)")
    for k, f in enumerate(fichas[:10], 1):
        print(f"  {k:>2}. {f['nombre'][:40]:<40} índice {fmt(f['indice']):>5}  cobertura {f['cobertura']*100:3.0f} %  {f['clase']}")


if __name__ == "__main__":
    main()
