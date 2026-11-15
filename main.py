#These codes were partially developed with Claude AI

#This code is used within the Python and QGIS console; 
# Purpose: Compute the threshold that defines the top 5% of pixels in the
# Mars Similarity Index raster (indice_marte.tif). Pixels at or above this
# value are considered the most Mars-like areas in Mexico and are used to
# build the candidate zones (Raster Calculator: "indice_marte@1" >= threshold).

import numpy as np
from osgeo import gdal

RUTA = r"C:\Users\sofia\OneDrive\Documentos\indice_marte.tif"

ds = gdal.Open(RUTA)
banda = ds.GetRasterBand(1)
a = banda.ReadAsArray().astype(float)

nodata = banda.GetNoDataValue()
if nodata is not None:
    a[a == nodata] = np.nan
a[(a < 0) | (a > 1)] = np.nan

umbral = np.nanpercentile(a, 95)
print("Umbral del 5% más alto:", round(umbral, 4))

#This code is used within the Python and QGIS console; 
#Purpose: its function was to provide the top 10 areas that are most similar to Mars 
#in Mexico based on data obtained from the analysis of the 4 variables 
#(precipitation, zero vegetation, dust, and temperature).

import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ARCHIVO = r"C:\Users\sofia\OneDrive\Documentos\resultados_zonas_excel.csv"

CARPETA = os.path.dirname(ARCHIVO)
print("Leyendo:", ARCHIVO)

df = pd.read_csv(ARCHIVO)
print("Columnas encontradas:", list(df.columns))

df = df.sort_values("ind_mean", ascending=False).reset_index(drop=True)
df["ranking"] = df.index + 1

columnas = ["ranking", "area_km2", "ind_mean", "ndvi_mean", "lluvia_mean", "rango_mean", "polvo_mean"]
columnas = [c for c in columnas if c in df.columns]
print(df[columnas].head(10))

df.to_csv(os.path.join(CARPETA, "ranking_zonas.csv"), index=False)

ax = df.head(10).plot.bar(x="ranking", y="ind_mean", legend=False,
                          title="Top 10 zonas: índice de similitud con Marte")
ax.set_xlabel("Ranking")
ax.set_ylabel("Índice (0 a 1)")
plt.tight_layout()
plt.savefig(os.path.join(CARPETA, "top10.png"), dpi=200)
plt.close()
print("Listo: ranking_zonas.csv y top10.png guardados en", CARPETA)