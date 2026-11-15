# Aztroteca
Aztroteca is a search engine for locations on Earth that are analogous to Martian environments.
# Mars Analog Sites in Mexico – NASA Space Apps Challenge 2026

We built a Mars Similarity Index for Mexico by combining satellite-derived criteria
and identified the top candidate zones as terrestrial Mars analogs.

## How to open
1. Download or clone this repository.
2. Open the file in QGIS (version 3.40.10)
3. Unzip all the folders from "some files" and "other files" and combine them into one with the QGIS file (the reason there are several folders is due to the size of the files used in QGIS).

## Data sources for the index
Format: | Criterion | Source | Product | Resolution |

| Aridity | NASA Giovanni | GPM_3IMERGM v07 (precipitation, 2015–2024) | 0.1° |
| Dust | NASA Giovanni | M2TMNXAER v5.12.4 (dust surface mass concentration) | 0.5° × 0.625° |
| Vegetation | NASA AppEEARS | MOD13A3.061, _1_km_monthly_NDVI (2024) | 1 km |
| Thermal range | NASA AppEEARS | MOD11A2.061, LST_Day_1km / LST_Night_1km (2024) | 1 km |

## Method
1. Data cleaning: time averages, clip to Mexico, resample to a common 1 km grid, water removal (NDVI < 0).
2. Normalization of each criterion to 0–1 (1 = more Mars-like).
3. Index = 0.30·aridity + 0.30·vegetation + 0.25·thermal range + 0.15·dust.
4. Candidate zones: top 5% of index values, polygonized, filtered (> 25 km²).
5. Zonal statistics and ranking (Python and QGIS).

## Folders
 –Python-code-used-in-QGIS-Albireo – Python code used in QGIS.
- QGIS_and_some_filestif – QGIS Program and some files used in QGIS (.tif).
- Some_filestif_p2 – Some files used in QGIS (.tif).
- Some_filestif_p3 – Some files used in QGIS (.tif).
- Some_filestif_p4 – Some files used in QGIS (.tif). 
- Some_filestif_p5 – Some files used in QGIS (.tif).
- other_files_qgis – Some files used in QGIS.
- Results – excel files, ranking, top 10, description of the zones.
- Images_evidence_QGIS – Some images of the process in QGIS.

## Limitations
Coarse dust resolution, one year of NDVI and temperature data, and Earth cannot reproduce
Mars gravity, atmospheric pressure or radiation.

## Team name
Albireo

##Team members
Sofía Castillo de Santiago.
Iniria Corichi Guevara.
Fátima Karuna Menéndez Cruz.
Yeleinie Stephanie Vazquez Novelo.
