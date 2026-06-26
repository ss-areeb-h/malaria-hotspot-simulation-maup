# PIN Code-based Malaria Hotspot Simulation and Detection

This repository contains R code for simulating malaria hotspot detection accuracy at different administrative levels in India. The code generates artificial hotspots, detects clusters using SaTScan, and evaluating hotspot detection performance across multiple administrative levels (District, TEHSIL/Subdistrict, Pincode).

The code is designed to support research on spatial epidemiology and hotspot detection for addressing and quantifying the modifiable area unit problem, with a focus on reproducible simulation and systematic performance assessment.

---

## Features

- Simulation of equal grid of baseline malaria cases within Indian states based on average case counts.
- Generation of artificial spatial hotspots using regular grids and random sampling.
- Aggregation of simulated cases to District, TEHSIL (subdistrict), and Pincode levels.
- Hotspot detection using SaTScan (purely spatial Poisson model).
- Performance evaluation using confusion matrix metrics, Jaccard index, and area coverage.
- Export of GeoPackage and CSV files for use in GIS and statistical analysis software.

---

## Repository Overview

```malaria-hotspot-simulation-map
.
├── README.md                     # Project description and usage
├── LICENSE                       # License for the code (add your chosen license)
├── pin_hotspot_simulation.R      # Main R script (your code)

```

## Data and Input Requirements

### Working Directory

The script assumes a working directory:

```r
PATH_WD <- "ENTER WORKING DIRECTORY HERE"
OUT_DIR <- "outputs2"
```

You should set `PATH_WD` to the folder where your shapefiles and `avg_malaria.csv` reside. The script will create `OUT_DIR` if it does not exist.

### Shapefiles

The following shapefiles are required inside the `Shapefiles/` subdirectory:

- `STATE_BOUNDARY.shp`  
  - Must contain a `STATE` field (state name).
- `DISTRICT_BOUNDARY.shp`  
  - Must contain `STATE` and `District` (or equivalent) identifiers.
- `SUBDISTRICT_BOUNDARY_fixed.shp`  
  - Must contain `STATE`, `District`, and `TEHSIL` (subdistrict) identifiers.
- `pincodes_official_fixed.shp`  
  - Must contain `STATE`, `Pincode`, and appropriate linking fields to districts/subdistricts.

All shapefiles are transformed to WGS84 (`EPSG:4326`) and invalid geometries are repaired before analysis.

### Malaria Cases Table

`avg_malaria.csv` is expected to contain at least:

- `State.UTs`: state or union territory name (matching `STATE` in the shapefile).
- `avg`: average malaria case count for that state.

The script joins this table to the `state_sf` layer to determine the number of simulated cases per state.

---

## Methods and Workflow

### 1. Baseline Case Simulation

For each state with non-zero average cases (`avg > 0` and `avg > 50`), the script:

1. Computes the state’s bounding box.
2. Generates a regular grid of points within the bounding box, or random points if case numbers are small.
3. Filters points to those that fall within the state polygon.
4. Adjusts the grid to ensure at least the required number of points (`n_cases`) is obtained.
5. Tags each point with the state name, creating the baseline case distribution.

States with insufficient data are skipped.

### 2. Artificial Hotspot Creation

Within each state:

1. The script transforms baseline points to an appropriate UTM projection based on the state centroid.
2. It repeatedly samples random hotspot centers and radii (up to `max_radius`) within the projected bounding box.
3. For each candidate hotspot:
   - Ensures limited overlap (≤ 10%) with previously created hotspots.
   - Identifies baseline points inside the hotspot buffer.
   - Adds additional simulated cases within the hotspot, scaled by `intensity_multiplier`.
4. All additional hotspot points are transformed back to WGS84.
5. The script constructs “true” hotspot polygons from the hotspot buffers for later evaluation.

### 3. Aggregation to Administrative Levels

Baseline and modified (baseline + hotspot) points are spatially joined to:

- District polygons
- TEHSIL (subdistrict) polygons
- Pincode polygons

For each level, case counts are aggregated per administrative unit, producing:

- Baseline counts (used as population in SaTScan)
- Modified counts (used as observed cases)

### 4. Hotspot Detection with SaTScan

For each state and each level (District, TEHSIL, Pincode):

1. The script writes three SaTScan input files:
   - `cases.cas`: observed cases with a dummy date.
   - `population.pop`: baseline counts (population at risk) with a dummy date.
   - `coordinates.geo`: geographic coordinates per administrative unit.
2. It configures SaTScan to run a purely spatial Poisson model:
   - High-rate clusters only.
   - No temporal dimension.
   - Monte Carlo simulations (999 replicates).
3. SaTScan is invoked via `rsatscan`, and spatial clusters (`shapeclust`) are read back into R.
4. Detected hotspot polygons are saved as GeoPackages per state and level.

**Note:** The script assumes SaTScan is installed at `C:/Program Files/SaTScan` and uses `SaTScanBatch64`. Adjust these paths if your installation differs.

### 5. Performance Evaluation

For each state and level:

1. True hotspots (`true_hotspot_polygons`) are compared with detected clusters (`shapeclust`).
2. The script computes:
   - Intersection areas and proportion of true hotspot area covered.
   - Jaccard index for overlap between true and detected clusters.
3. At the administrative unit level, it flags each unit as:
   - `hot_ref`: intersects any true hotspot.
   - `hot_test`: intersects any detected cluster.
4. Confusion matrix metrics are computed:
   - True Positives (TP), False Negatives (FN), False Positives (FP), True Negatives (TN).
   - Sensitivity, specificity, positive predictive value (PPV), negative predictive value (NPV).
   - Mean proportion of hotspot area covered and mean Jaccard index.
5. These metrics are stored per state, level, and simulation.

---

## Outputs

After `n_simulations` runs:

- **Per simulation and per state**:
  - `STATE_baseline_pts.gpkg`: baseline simulated case points.
  - `STATE_hotspot_pts.gpkg`: hotspot-only points.
  - `STATE_hotspot_poly.gpkg`: true hotspot polygons.
  - `STATE_hotspot_district.gpkg`: detected clusters at district level (if any).
  - `STATE_hotspot_tehsil.gpkg`: detected clusters at TEHSIL level (if any).
  - `STATE_hotspot_pincode.gpkg`: detected clusters at Pincode level (if any).

- **Global CSV outputs** (in `PATH_WD`):
  - `simulation_results_detailed.csv`  
    - Accuracy metrics for each state, level, and simulation.
  - `simulation_results_summary.csv`  
    - Average metrics per state and level across all simulations.

These outputs can be used for statistical analysis, visualization in GIS software, and reporting.

---

## Requirements

- **R**: Version 4.0 or later (tested on recent R versions).
- **R packages**:
  - `sf`
  - `sp`
  - `terra`
  - `dplyr`
  - `spatstat`
  - `rsatscan`
  - `spdep`
  - `purrr`
  - `exactextractr`
- **SaTScan**:
  - SaTScan installed locally.
  - Accessible batch executable (`SaTScanBatch64`) and parameter files.
  - Adjust `sslocation` and `ssbatchfilename` in the script if your installation differs.

Install the required R packages with:

```r
install.packages(c(
  "sf", "sp", "terra", "dplyr", "spatstat",
  "spdep", "purrr", "exactextractr"
))
# rsatscan may require installation from a specific source (CRAN or other)
install.packages("rsatscan")
```

---

## How to Run

1. Clone or download this repository from GitHub.
2. Open `pin_hotspot_simulation.R` in R or RStudio.
3. Adjust `PATH_WD`, `OUT_DIR`, and any file paths to match your local directory structure.
4. Ensure all required shapefiles and `avg_malaria.csv` are present.
5. Confirm SaTScan is installed and paths in the script are correct.
6. Run the script. It will:
   - Perform up to `n_simulations` simulations.
   - Create subfolders under `OUT_DIR` for each simulation.
   - Save per-state GeoPackages and global CSV results.

Depending on the number of simulations and size of your spatial data, runtime may vary from several minutes to longer.

---

## Reproducibility and Extensions

- You can change key parameters at the top of the script:
  - `n_simulations`
  - `n_hotspots`
  - `max_radius`
  - `intensity_multiplier`
- You can adapt the code to:
  - Different diseases or case distributions.
  - Alternative administrative hierarchies or boundary datasets.
  - Different SaTScan settings (e.g., adding temporal dimension, changing model type).

For reproducible research, consider:
- Archiving the exact version of this script alongside input data.
- Using Git tags or releases to denote versions tied to specific analyses.
- Publishing the dataset and code in a Dataverse installation to obtain a DOI.

---

## Citation

If you use this code in academic work, please cite it as:

> Hussain SSA; Yadav CP, Sharma A. Year. *PIN Code-based malaria hotspot simulation and detection*. GitHub repository: `<GitHub URL>`. Dataverse dataset: `<DOI>`.

---

## License

> This project is licensed under the MIT License – see the `LICENSE` file for details.
