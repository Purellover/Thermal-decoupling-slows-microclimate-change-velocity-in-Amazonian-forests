# Analysis and figure generation scripts

This repository contains simplified Python scripts used to train the microclimate models, calculate climate change velocity, identify microclimate refugia and generate the main figures for the manuscript *Thermal decoupling slows microclimate change velocity in Amazonian forests*.

The scripts assume that the required public climate, vegetation, topographic and intermediate raster/vector files have already been prepared following the processing workflow described in the manuscript. No local machine-specific paths are used. Paths are defined relative to the repository root.

## Repository structure

```text
repository/
├── data/
│   ├── model_input/
│   │   └── Merged_Data.csv
│   ├── shapefiles/
│   │   ├── amazon_bio.shp
│   │   └── hotpoint.shp
│   ├── fig1/
│   ├── fig2/
│   ├── fig3/
│   ├── fig4/
│   ├── historical_baseline/
│   ├── terrain_structure_data/
│   ├── prediction_annual_mean_tmax/
│   └── prediction_annual_mean_tmin/
├── scripts/
├── outputs/
├── requirements.txt
└── README.md
```

## Scripts

`scripts/train_microclimate_xgboost_models.py` trains the XGBoost models used to predict understory mean, maximum and minimum air temperature. The script includes feature construction, random cross-validation, within-site validation, leave-one-site-out spatial block validation, feature importance analysis and partial dependence plots.

`scripts/calculate_climate_velocity.py` calculates horizontal climate change velocity from temporal warming trends and historical spatial temperature gradients. The script estimates velocity magnitude and displacement direction for microclimate and macroclimate layers and supports sensitivity testing across spatial gradient resolutions.

`scripts/make_figure1_warming_decoupling.py` generates warming-rate maps and ΔWR summaries across environmental gradients.

`scripts/make_figure2_microclimate_regimes.py` generates thermal regime maps, Sankey transitions and regime thermal characteristics. Regimes are fitted on the historical baseline and future periods are assigned to the same historical cluster centroids.

`scripts/make_figure3_climate_velocity.py` generates microclimate and macroclimate velocity maps with hotspot rose diagrams.

`scripts/make_figure4_refugia_protection.py` generates refugia maps, climate residence time, area and patch-count summaries, and a supplementary protected-area gap figure.

Merged_Data.csv is the prepared model training table containing matched monthly microclimate observations and predictor variables. The original microclimate logger data are not distributed in this repository because of field-site data sharing agreements, but they are available from the corresponding author upon reasonable request.

## System requirements

The scripts were developed and tested with Python 3.10 on a Windows desktop environment. They should also run on Linux or macOS provided that the required geospatial Python libraries are correctly installed.

The main Python dependencies are listed in `requirements.txt` and include `numpy`, `pandas`, `matplotlib`, `rasterio`, `geopandas`, `shapely`, `scikit-learn`, `xgboost`, `scipy`, `plotly`, `Pillow` and `PyMuPDF`.

No non-standard hardware is required. A desktop computer with at least 16 GB RAM is recommended for processing large raster files. Some basin-wide raster operations may benefit from higher memory availability.

## Installation

Create a Python environment and install the dependencies:

```bash
conda create -n microclimate_velocity python=3.10
conda activate microclimate_velocity
pip install -r requirements.txt
```

Typical installation time on a standard desktop computer is approximately 5–15 minutes, depending on internet speed and whether geospatial dependencies are already available.

Some Sankey or PDF-image embedding steps may require additional system dependencies depending on the operating system, such as Kaleido for Plotly static export or PyMuPDF for reading PDF panels.

## Input data

The scripts require prepared tabular, raster and vector files generated during the analysis workflow. These files should be placed in the `data/` directory following the structure shown above.

Large raster datasets and restricted field-derived intermediate products are not included directly in this repository because of file size and data-sharing constraints. Publicly available input datasets are described in the manuscript Data Availability section. The microclimate logger data are available from the corresponding author upon reasonable request, subject to field-site data sharing agreements.

## Running the scripts

Run the scripts from the repository root. A typical workflow is:

```bash
python scripts/train_microclimate_xgboost_models.py
python scripts/calculate_climate_velocity.py
python scripts/make_figure1_warming_decoupling.py
python scripts/make_figure2_microclimate_regimes.py
python scripts/make_figure3_climate_velocity.py
python scripts/make_figure4_refugia_protection.py
```

The figure scripts can also be run independently if the required intermediate raster and vector files are already available.

## Expected outputs

The scripts write model outputs, velocity rasters, figure files and summary tables to the `outputs/` directory. Expected outputs include:

```text
outputs/
├── model_outputs/
│   ├── xgboost_Tmean_model.pkl
│   ├── xgboost_Tmax_model.pkl
│   ├── xgboost_Tmin_model.pkl
│   ├── validation_metrics.csv
│   ├── loso_site_metrics.csv
│   └── feature_importance_tables/
├── velocity_outputs/
│   ├── climate_velocity_ssp245_10km.tif
│   ├── climate_velocity_ssp585_10km.tif
│   ├── velocity_direction_ssp245_10km.tif
│   └── velocity_direction_ssp585_10km.tif
├── Fig1_SSP245.pdf
├── Fig1_SSP585.pdf
├── Fig2_SSP245.pdf
├── Fig2_SSP585.pdf
├── Fig3_SSP245.pdf
├── Fig3_SSP585.pdf
├── Fig4_refugia_analysis.pdf
├── Fig4_supplementary_gap.pdf
└── related CSV or shapefile summary outputs
```

Output filenames may vary slightly depending on the script settings.

## Expected runtime

Runtime depends on raster size and local hardware. On a standard desktop computer, model training typically requires several minutes to tens of minutes, depending on the size of the training table and the hyperparameter search settings. Climate velocity calculation and figure generation can take from several minutes to less than one hour per script. Figure 2 and Figure 4 may take longer because they include clustering, Sankey diagram generation, protected-area overlay or raster-vector operations.

## Demo data

A small demo dataset is not included in this repository. The scripts are intended to reproduce the manuscript analyses using the full prepared tabular, raster and vector datasets generated by the analysis workflow. Users can test the scripts by providing reduced versions of the required input files using the same file names and directory structure.

## Reproduction instructions

To reproduce the main outputs, place the required input and intermediate files in the `data/` directory, install the dependencies listed above, and run the scripts from the repository root. The resulting model diagnostics, climate velocity rasters, figures and summary outputs will be saved to the `outputs/` directory.
