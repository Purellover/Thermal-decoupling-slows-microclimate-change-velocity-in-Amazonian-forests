# Figure generation scripts

This repository contains simplified Python scripts used to generate the main figures for the manuscript
*Thermal decoupling slows microclimate change velocity in Amazonian forests*.

The scripts assume that the required intermediate raster and vector files have already been generated.
No local machine-specific paths are used. Paths are defined relative to the repository root:

```text
repository/
├── data/
│   ├── shapefiles/
│   │   ├── amazon_bio.shp
│   │   └── hotpoint.shp
│   ├── fig1/
│   │   ├── Trend_Micro_Tmean_ssp245_1983_2100.tif
│   │   ├── Trend_Macro_Tmean_Masked_ssp245_1983_2100.tif
│   │   ├── Delta_WarmingRate_ssp245_1983_2100.tif
│   │   └── ...
│   ├── fig3/
│   │   ├── ssp245/10km/climate_velocity_ssp245_10km.tif
│   │   ├── macro/ssp245/10km/climate_velocity_ssp245_10km.tif
│   │   └── rose_hotspot/
│   ├── fig4/
│   │   ├── ssp245/core_ssp245.shp
│   │   ├── ssp585/core_ssp585.shp
│   │   └── protected_areas/
│   ├── historical_baseline/
│   ├── terrain_structure_data/
│   ├── prediction_annual_mean_tmax/
│   └── prediction_annual_mean_tmin/
├── scripts/
└── outputs/
```

## Scripts

- `scripts/make_figure1_warming_decoupling.py` generates warming-rate maps and ΔWR summaries across environmental gradients.
- `scripts/make_figure2_microclimate_regimes.py` generates thermal regime maps, Sankey transitions and regime thermal characteristics. Regimes are fitted on the historical baseline and future periods are assigned to the same historical cluster centroids.
- `scripts/make_figure3_climate_velocity.py` generates microclimate and macroclimate velocity maps with hotspot rose diagrams.
- `scripts/make_figure4_refugia_protection.py` generates refugia maps, climate residence time, area and patch-count summaries, and a supplementary protected-area gap figure.

## Installation

```bash
pip install -r requirements.txt
```

Some Sankey or PDF-image embedding steps may require additional system dependencies depending on the operating system, such as Kaleido for Plotly static export or PyMuPDF for reading PDF panels.

## Running the scripts

Run from the repository root:

```bash
python scripts/make_figure1_warming_decoupling.py
python scripts/make_figure2_microclimate_regimes.py
python scripts/make_figure3_climate_velocity.py
python scripts/make_figure4_refugia_protection.py
```

Outputs are written to the `outputs/` directory.
