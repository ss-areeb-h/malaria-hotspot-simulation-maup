# Load required libraries
library(sf)           # Spatial data handling
library(sp)
library(terra)        # Raster operations
library(dplyr)        # Data manipulation
library(spatstat)     # Point pattern analysis
library(rsatscan)     # SaTScan integration
library(spdep)        # Spatial autocorrelation
library(purrr)        # Functional programming
library(exactextractr)


# 2. User-defined inputs: paths and field names
PATH_WD <- "C:/Users/ss_ar/OneDrive/Work/PIN code/Simulation/Equal Grids - States"
PATH_STATE_SHP <- "Shapefiles/STATE_BOUNDARY.shp"
PATH_DIST_SHP <- "Shapefiles/DISTRICT_BOUNDARY.shp"             # must have: district_id
PATH_SUBDIST_SHP <- "Shapefiles/SUBDISTRICT_BOUNDARY_fixed.shp"       # must have: subdistrict_id, district_id
PATH_PIN_SHP <- "Shapefiles/pincodes_official_fixed.shp"               # must have: pincode, subdistrict_id, district_id
OUT_DIR <- "outputs"

## Working and output directory
setwd(PATH_WD)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## Load Files
state_sf <- st_read(PATH_STATE_SHP, quiet = TRUE)
dist_sf <- st_read(PATH_DIST_SHP, quiet = TRUE)
subdist_sf <- st_read(PATH_SUBDIST_SHP, quiet = TRUE)
pin_sf <- st_read(PATH_PIN_SHP, quiet = TRUE)
mal_cases <- read.csv("avg_malaria.csv")

### subdistrict adjustments
subdist_sf <- subdist_sf %>%
  mutate(
    TEHSIL = if_else(
      is.na(TEHSIL),          # where TEHSIL is missing
      paste0(District),  # District value + "2" at the end
      TEHSIL                  # otherwise keep existing TEHSIL
    )
  )

subdist_sf <- subdist_sf %>%
  group_by(TEHSIL) %>%
  mutate(
    TEHSIL = if(n() > 1) paste0(TEHSIL, 1:n()) else TEHSIL
  ) %>%
  ungroup()
subdist_sf <- subdist_sf[!st_is_empty(subdist_sf), ]

## Align Coordinate Reference Systems
TARGET_CRS <- 4326
state_sf <- st_transform(state_sf, TARGET_CRS)
dist_sf <- st_transform(dist_sf, TARGET_CRS)
subdist_sf <- st_transform(subdist_sf, TARGET_CRS)
pin_sf <- st_transform(pin_sf, TARGET_CRS)

state_sf <- st_make_valid(state_sf)
dist_sf <- st_make_valid(dist_sf)
subdist_sf <- st_make_valid(subdist_sf)
pin_sf <- st_make_valid(pin_sf)

# Initialize parameters for simulation
n_simulations <- 10
n_hotspots <- 5
max_radius <- 50000
intensity_multiplier <- 10

# Joining Malaria Cases to the States shapefile
state_sf_mal <- state_sf %>% left_join(mal_cases, by=c("STATE" = "State.UTs"))

# Initialize results storage
all_sim_results <- list()

# DISABLE S2 GLOBALLY
sf_use_s2(FALSE)

# START SIMULATION LOOP
for(sim in 1:n_simulations) {
  sim=1
  cat("Running simulation", sim, "of", n_simulations, "\n")
  
  # Create simulation-specific output directory
  sim_dir <- file.path(OUT_DIR, paste0("simulation_", sim))
  dir.create(sim_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Remove states with 0 or NA cases
  valid_states <- state_sf_mal[!is.na(state_sf_mal$avg) & state_sf_mal$avg > 0, ]
  cat("Processing", nrow(valid_states), "states with cases\n")
  excluded_states <- data.frame("State"=character(), "No. of Cases" = integer(), stringsAsFactors = FALSE)
  
  sim_accuracy_metrics <- list()
  
  # STATE LOOP
  for(state_idx in 20:nrow(valid_states)) {
    state_idx=7
    state <- valid_states[state_idx, ]
    state_name <- state$STATE
    n_cases <- round(state$avg, 0)
    
    cat("Processing state", state_name, "with", n_cases, "cases\n")
    
    
    if(n_cases <= 50) {
      excluded_states <- rbind(excluded_states,
                               data.frame(State = state_name,
                                          `No. of Cases` = n_cases,
                                          stringsAsFactors = FALSE))
      next
    } else {
    # Get state-specific administrative boundaries
    districts <- dist_sf[dist_sf$STATE == state_name, ]
    subdistricts <- subdist_sf[subdist_sf$STATE == state_name, ]
    pincodes <- pin_sf[pin_sf$STATE == state_name, ]
    
    # Check if we have administrative data
    if(nrow(districts) == 0 || nrow(subdistricts) == 0 || nrow(pincodes) == 0) {
      cat("  Skipping state", state_name, "- missing administrative boundaries\n")
      next
    }
    
    admin_sf <- list(District = districts,
                     TEHSIL = subdistricts,
                     Pincode = pincodes)
    
    ## SECTION 1: GENERATE EQUIDISTANT POINTS WITHIN STATES BASED ON NO. OF CASES ##  
    
    baseline_points <- NULL
    tryCatch({
      # Create a regular grid of points within the state boundary
      state_bbox <- st_bbox(state)
      
      # Calculate approximate spacing for equidistant points
      bbox_width <- state_bbox["xmax"] - state_bbox["xmin"]
      bbox_height <- state_bbox["ymax"] - state_bbox["ymin"]
      area <- bbox_width * bbox_height
      
      # If we have very few points, use random sampling instead of grid
      if(n_cases < 10) {
        points_sf <- st_sample(state, size = n_cases, type = "random")
        points_sf <- st_as_sf(points_sf)
        colnames(points_sf)[1] <- "geometry" 
      } else {
        # Calculate optimal grid dimensions
        aspect_ratio <- bbox_width / bbox_height
        n_cols <- round(sqrt(n_cases * aspect_ratio))
        n_rows <- round(n_cases / n_cols)
        
        # Ensure we have at least the required number of points
        while(n_cols * n_rows < n_cases) {
          n_cols <- n_cols + 1
          n_rows <- round(n_cases / n_cols)
        }
        
        # Create grid points
        x_seq <- seq(state_bbox["xmin"], state_bbox["xmax"], length.out = n_cols + 2)[2:(n_cols + 1)]
        y_seq <- seq(state_bbox["ymin"], state_bbox["ymax"], length.out = n_rows + 2)[2:(n_rows + 1)]
        
        grid_points <- expand.grid(x = x_seq, y = y_seq)
        grid_sf <- st_as_sf(grid_points, coords = c("x", "y"), crs = st_crs(state))
        
        # Keep only points that fall within the state boundary
        points_within <- st_intersects(grid_sf, state)
        valid_indices <- which(sapply(points_within, length) > 0)
        
        if(length(valid_indices) >= n_cases) {
          # Sample the required number of points from valid grid points
          sampled_indices <- sample(valid_indices, size = n_cases)
          points_sf <- grid_sf[sampled_indices, ]
        } else {
          # If not enough grid points, supplement with random points
          points_sf <- grid_sf[valid_indices, ]
          additional_needed <- n_cases - nrow(points_sf)
          if(additional_needed > 0) {
            additional_points <- st_sample(state, size = additional_needed, type = "random")
            additional_sf <- st_as_sf(additional_points)
            names(additional_sf)[names(additional_sf) == "x"] <- "geometry"
            additional_sf <- st_set_geometry(additional_sf, "geometry")
            points_sf <- rbind(points_sf, additional_sf)
          }
        }
      }
      
      # Add state name to points
      points_sf$state <- state_name
      baseline_points <- points_sf[, c("state", "geometry")]
      
      cat("  Successfully generated", n_cases, "points for state", state_name, "\n")
      
    }, error = function(e) {
      cat("  ERROR in State", state_name, ":", e$message, "\n")
      # Fallback: use simple random sampling
      try({
        points_sf <- st_sample(state, size = n_cases, type = "random")
        points_sf <- st_as_sf(points_sf)
        points_sf$state <- state_name
        baseline_points <- points_sf[, c("state", "geometry")]
        cat("  Used fallback random sampling for State", state_name, "\n")
      })
    })
    
    if(is.null(baseline_points)) next
    
    ## SECTION 2: CREATE ARTIFICIAL HOTSPOTS ##
    
    modified_points <- baseline_points
    hotspot_info <- list()
    
    # Get appropriate UTM zone for this state
    centroid <- st_centroid(state)
    coords <- st_coordinates(centroid)
    lon <- coords[1, "X"]
    lat <- coords[1, "Y"]
    utm_zone <- floor((lon + 180) / 6) + 1
    if (lat >= 0) {
      utm_crs <- 32600 + utm_zone
    } else {
      utm_crs <- 32700 + utm_zone
    }
    
    # Transform to UTM CRS
    baseline_points_proj <- st_transform(baseline_points, utm_crs)
    modified_points_proj <- st_transform(modified_points, utm_crs)
    
    # Get bounding box of study area in projected CRS
    bbox <- st_bbox(baseline_points_proj)
    
    hotspots_created <- 0
    attempts <- 0
    max_attempts <- 200
    
    while(hotspots_created < n_hotspots && attempts < max_attempts) {
      attempts <- attempts + 1
      
      # Random hotspot center in projected CRS
      center_x <- runif(1, bbox["xmin"], bbox["xmax"])
      center_y <- runif(1, bbox["ymin"], bbox["ymax"])
      center <- st_point(c(center_x, center_y))
      center_sf <- st_sfc(center, crs = utm_crs)
      
      # Random radius (10% to 100% of max_radius)
      radius <- runif(1, max_radius * 0.1, max_radius)
      
      # Create circular buffer around center in projected CRS
      hotspot_buffer_proj <- st_buffer(center_sf, dist = radius)
      
      # Check for overlap with existing hotspots (max 10% overlap allowed)
      overlap_too_high <- FALSE
      if(length(hotspot_info) > 0) {
        for(existing_hotspot in hotspot_info) {
          existing_poly <- existing_hotspot$geometry_proj
          intersection <- st_intersection(hotspot_buffer_proj, existing_poly)
          if(length(intersection) > 0 && !inherits(intersection, "sfc_GEOMETRY")) {
            intersection_area <- as.numeric(st_area(intersection))
            new_area <- as.numeric(st_area(hotspot_buffer_proj))
            overlap_ratio <- intersection_area / new_area
            if(overlap_ratio > 0.1) {
              overlap_too_high <- TRUE
              break
            }
          }
        }
      }
      
      if(overlap_too_high) {
        next
      }
      
      # Find points within hotspot area
      points_in_hotspot <- st_intersects(hotspot_buffer_proj, baseline_points_proj)[[1]]
      
      if(length(points_in_hotspot) > 10) {
        hotspots_created <- hotspots_created + 1
        
        # Calculate additional cases (intensity multiplier)
        n_additional_cases <- round(length(points_in_hotspot) *
                                      runif(1, intensity_multiplier * 0.5, intensity_multiplier))
        
        # Generate additional points within hotspot in projected CRS
        additional_points_proj <- st_sample(hotspot_buffer_proj, size = n_additional_cases, type = "random")
        
        # Convert to sf and ensure proper column structure
        additional_df_proj <- st_as_sf(additional_points_proj)
        st_geometry(additional_df_proj) <- "geometry"
        additional_df_proj$state <- state_name
        
        # Transform back to geographic CRS
        additional_df <- st_transform(additional_df_proj, 4326)
        
        # Add to modified points
        modified_points_proj <- rbind(
          modified_points_proj[, c("state", "geometry")],
          additional_df_proj[, c("state", "geometry")]
        )
        
        # Store hotspot metadata
        hotspot_buffer_geo <- st_transform(hotspot_buffer_proj, 4326)
        
        hotspot_info[[hotspots_created]] <- list(
          center = center,
          radius = radius,
          original_cases = length(points_in_hotspot),
          additional_cases = n_additional_cases,
          geometry_proj = hotspot_buffer_proj,
          geometry = hotspot_buffer_geo
        )
      }
    }
    
    # Convert all points back to geographic CRS
    modified_points <- st_transform(modified_points_proj, 4326)
    
    # Extract hotspot points
    baseline_coords <- st_coordinates(baseline_points)
    modified_coords <- st_coordinates(modified_points)
    is_hotspot <- !paste(modified_coords[,1], modified_coords[,2]) %in% 
      paste(baseline_coords[,1], baseline_coords[,2])
    hotspot_points <- modified_points[is_hotspot, ]
    
    # Method 1: Using st_within to filter points inside the state
    points_within_state <- st_within(hotspot_points, state)
    valid_points_indices <- which(sapply(points_within_state, length) > 0)
    hotspot_points <- hotspot_points[valid_points_indices, ]
    
    # Also ensure modified_points are within state boundary
    points_within_state_mod <- st_within(modified_points, state)
    valid_points_indices_mod <- which(sapply(points_within_state_mod, length) > 0)
    modified_points <- modified_points[valid_points_indices_mod, ]
    
    # Create true hotspot polygons
    if(length(hotspot_info) > 0) {
      true_hotspot_polygons <- do.call(rbind, lapply(hotspot_info, function(hs) {
        st_sf(geometry = st_sfc(hs$geometry), crs = 4326)
      }))
    } else {
      true_hotspot_polygons <- st_sf(geometry = st_sfc(), crs = 4326)
    }
    
    cat("Total no. of cases in all hotspots:", nrow(hotspot_points), "\n")
    cat("Total points after adding hotspots:", nrow(modified_points), "\n")
    
    if (nrow(hotspot_points) == 0) {
      excluded_states <- rbind(
        excluded_states,
        data.frame(
          State = state_name,
          `No. of Cases` = n_cases,
          stringsAsFactors = FALSE
        )
      )
      next  # skip to next state
    }
    
    
    ## SECTION 3: AGGREGATE TO ADMINISTRATIVE LEVELS ##
    
    # Aggregating baseline points to administrative levels for calibration
    district_case_points_baseline <- baseline_points %>% st_join(districts[, c("District")], join = st_within)
    district_counts_baseline <- district_case_points_baseline %>%
      st_drop_geometry() %>%
      group_by(District) %>%
      summarise(cases = n()) %>%
      left_join(districts, ., by = "District")
    
    subdistrict_case_points_baseline <- baseline_points %>% st_join(subdistricts[, c("TEHSIL")], join = st_within)
    subdistrict_counts_baseline <- subdistrict_case_points_baseline %>%
      st_drop_geometry() %>%
      group_by(TEHSIL) %>%
      summarise(cases = n()) %>%
      left_join(subdistricts, ., by = "TEHSIL")
    
    pincode_case_points_baseline <- baseline_points %>% st_join(pincodes[, c("Pincode")], join = st_within)
    pincode_counts_baseline <- pincode_case_points_baseline %>%
      st_drop_geometry() %>%
      group_by(Pincode) %>%
      summarise(cases = n()) %>%
      left_join(pincodes, ., by = "Pincode")
    
    baseline_aggregated <- list(District = district_counts_baseline,
                                TEHSIL = subdistrict_counts_baseline,
                                Pincode = pincode_counts_baseline)
    
    # Aggregating modified points to administrative levels for hotspot detection
    district_case_points <- modified_points %>% st_join(districts[, c("District")], join = st_within)
    district_counts <- district_case_points %>%
      st_drop_geometry() %>%
      group_by(District) %>%
      summarise(cases = n()) %>%
      left_join(districts, ., by = "District")
    
    subdistrict_case_points <- modified_points %>% st_join(subdistricts[, c("TEHSIL")], join = st_within)
    subdistrict_counts <- subdistrict_case_points %>%
      st_drop_geometry() %>%
      group_by(TEHSIL) %>%
      summarise(cases = n()) %>%
      left_join(subdistricts, ., by = "TEHSIL")
    
    pincode_case_points <- modified_points %>% st_join(pincodes[, c("Pincode")], join = st_within)
    pincode_counts <- pincode_case_points %>%
      st_drop_geometry() %>%
      group_by(Pincode) %>%
      summarise(cases = n()) %>%
      left_join(pincodes, ., by = "Pincode")
    
    aggregated_data <- list(
      District = district_counts,
      TEHSIL = subdistrict_counts,
      Pincode = pincode_counts
    )
    
    # Create state simulation directory
    state_sim_dir <- file.path(sim_dir, state_name)
    dir.create(state_sim_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Save point files
    # Save as GeoPackage (a single .gpkg file) instead of Shapefile
    st_write(baseline_points, 
             file.path(state_sim_dir, paste0(state_name, "_baseline_pts.gpkg")),
             delete_dsn = TRUE, quiet = TRUE)
    
    st_write(hotspot_points, 
             file.path(state_sim_dir, paste0(state_name, "_hotspot_pts.gpkg")),
             delete_dsn = TRUE, quiet = TRUE)
    
    st_write(true_hotspot_polygons,
             file.path(state_sim_dir, paste0(state_name, "_hotspot_poly.gpkg")),
             delete_dsn = TRUE, quiet = TRUE)
    
    ## SECTION 4: HOTSPOT DETECTION AND PERFORMANCE EVALUATION BY LEVEL ##
    
    # LEVEL LOOP
    for(level in c("District", "TEHSIL", "Pincode")) {    ####"TEHSIL",
      cat("  Processing level:", level, "for state", state_name, "\n")
      
      # Prepare data for SaTScan
      baseline_data <- baseline_aggregated[[level]] %>%
        st_drop_geometry() %>%
        select(level, "cases")
      baseline_data$cases <- as.double(baseline_data$cases)
      baseline_data[is.na(baseline_data)] <- 0.1
      
      pop_data <- baseline_data %>% rename(pop = cases)
      pop_data[,"date"] <- "2000/01/01" 
      pop_data <- pop_data %>% relocate(date, .before = 2)
      
      obs_data <- aggregated_data[[level]] %>%
        st_drop_geometry() %>%
        select(level, "cases") 
      obs_data[is.na(obs_data)] <- 0
      obs_data[,"date"] <- "2000/01/01"
      
      coords <- st_point_on_surface(aggregated_data[[level]]) %>% st_coordinates() %>% as.data.frame()
      coords[,level] <- obs_data[[level]]
      coords <- coords %>% relocate(last_col(), .before = 1) %>% relocate(last_col(), .before = 2)
      
      # Write files for SaTScan
      write.table(obs_data, "cases.cas", sep = "\t", row.names = FALSE, col.names = FALSE)
      write.table(coords, "coordinates.geo", sep = "\t", row.names = FALSE, col.names = FALSE)
      write.table(pop_data, "population.pop", sep = "\t", row.names = FALSE, col.names = FALSE)
      
      # Run SaTScan
      satscan_results <- NULL
      tryCatch({
        
        library(rsatscan)
        
        # Set SaTScan parameters
        
        invisible(ss.options(reset = TRUE))
        ss.options(list(CaseFile="cases.cas",
                        PrecisionCaseTimes=1,
                        PopulationFile="population.pop", 
                        CoordinatesFile="coordinates.geo",
                        CoordinatesType=1,
                        AnalysisType=1,  # Purely spatial
                        ModelType=0,     # Poisson
                        ScanAreas=1,     # High rates only
                        TimeAggregationUnits=0,  # No time aggregation
                        MaxSpatialSizeInPopulationAtRisk=100,  # Max 50% of population
                        UseDistanceFromCenterOption="y",
                        MaxSpatialSizeInDistanceFromCenter=100,
                        MinimumCasesInHighRateClusters=20,
                        MonteCarloReps=999,
                        OutputGoogleEarthKML="n",
                        OutputShapefiles="y",
                        ResultsFile="results.txt"
        ))
        
        write.ss.prm(PATH_WD, "parameters")
        
        satscan_results <- satscan(prmlocation = PATH_WD, 
                                   prmfilename = "parameters",
                                   sslocation = "C:/Program Files/SaTScan",
                                   ssbatchfilename = "SaTScanBatch64",
                                   cleanup = TRUE,
                                   verbose = FALSE)
        
      }, error = function(e) {
        cat("    SaTScan error:", e$message, "\n")
        satscan_results <- NULL
      })
      
      # Clean up temporary files
      if(file.exists("cases.cas")) file.remove("cases.cas")
      if(file.exists("coordinates.geo")) file.remove("coordinates.geo")
      if(file.exists("population.pop")) file.remove("population.pop")
      if(file.exists("parameters.prm")) file.remove("parameters.prm")
      
      # Save detected hotspots
      if(!is.null(satscan_results) && !is.null(satscan_results$shapeclust)) {
        st_write(satscan_results$shapeclust,
                 file.path(state_sim_dir, paste0(state_name, "_hotspot_", tolower(level), ".gpkg")),
                 delete_dsn = TRUE, quiet = TRUE)
      }
      
      # Evaluate hotspots
      if(is.null(satscan_results) || is.null(satscan_results$shapeclust)) {
        # No hotspots detected
        accuracy_metrics <- data.frame(
          TP = 0, FN = 0, FP = 0, TN = 0,
          sensitivity = 0, specificity = 0, PPV = 0, NPV = 0,
          proportion_area_covered = 0, jaccard_avg = 0
        )
      } else {
        shp1 <- true_hotspot_polygons
        shp2 <- satscan_results$shapeclust
        
        # Create cluster IDs if they don't exist
        if (!"cluster_id" %in% names(shp1)) shp1$cluster_id <- seq_len(nrow(shp1))
        if (!"CLUSTER" %in% names(shp2)) shp2$CLUSTER <- seq_len(nrow(shp2))
        
        # Calculate intersection areas
        inter12 <- st_intersection(
          shp1 %>% select(cluster_id, geometry),
          shp2 %>% select(cluster_id2 = CLUSTER, geometry)
        )
        
        if(nrow(inter12) > 0) {
          inter12 <- inter12 %>%
            mutate(area_int = as.numeric(st_area(geometry)))
        } else {
          inter12 <- inter12 %>%
            mutate(area_int = numeric(0))
        }
        
        # Total area per reference cluster
        ref_area <- shp1 %>%
          mutate(area_ref = as.numeric(st_area(geometry))) %>%
          st_drop_geometry() %>%
          select(cluster_id, area_ref)
        
        # Summarise overlap per reference cluster
        if(nrow(inter12) > 0) {
          overlap_ref <- inter12 %>%
            st_drop_geometry() %>%
            group_by(cluster_id) %>%
            summarise(area_int_sum = sum(area_int, na.rm = TRUE), .groups = "drop") %>%
            right_join(ref_area, by = "cluster_id") %>%
            mutate(
              area_int_sum = if_else(is.na(area_int_sum), 0, area_int_sum),
              prop_ref_covered = area_int_sum / area_ref
            )
        } else {
          overlap_ref <- ref_area %>%
            mutate(area_int_sum = 0, prop_ref_covered = 0)
        }
        
        # Jaccard index calculation
        jaccard_ref <- map_dfr(shp1$cluster_id, function(cid) {
          ref_poly <- shp1 %>% filter(cluster_id == 1)
          test_over <- st_intersection(ref_poly, shp2)
          if (nrow(test_over) == 0) {
            return(tibble(cluster_id = cid, jaccard = 0))
          }
          test_union <- shp2[shp2$CLUSTER %in% unique(test_over$cluster_id), ] %>%
            st_union()
          test_union <- st_sf(geometry=test_union)
          inter_area <- as.numeric(st_area(st_intersection(ref_poly, test_union)))
          union_area <- as.numeric(st_area(st_union(ref_poly, test_union)))
          tibble(cluster_id = cid, jaccard = ifelse(union_area > 0, inter_area / union_area, 0))
        })
        
        # Merge overlap metrics
        cluster_metrics <- overlap_ref %>%
          left_join(jaccard_ref, by = "cluster_id")
        
        # Flag admin units as hotspot or not
        admin_sf[[level]]$hot_ref <- lengths(st_intersects(admin_sf[[level]], shp1)) > 0
        admin_sf[[level]]$hot_test <- lengths(st_intersects(admin_sf[[level]], shp2)) > 0
        
        # Confusion matrix components
        TP <- sum(admin_sf[[level]]$hot_ref & admin_sf[[level]]$hot_test)
        FN <- sum(admin_sf[[level]]$hot_ref & !admin_sf[[level]]$hot_test)
        FP <- sum(!admin_sf[[level]]$hot_ref & admin_sf[[level]]$hot_test)
        TN <- sum(!admin_sf[[level]]$hot_ref & !admin_sf[[level]]$hot_test)
        
        # Metrics
        sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
        specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), 0)
        PPV <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
        NPV <- ifelse((TN + FN) > 0, TN / (TN + FN), 0)
        proportion_area_covered <- mean(cluster_metrics$prop_ref_covered, na.rm = TRUE)
        jaccard_avg <- mean(cluster_metrics$jaccard, na.rm = TRUE)
        
        accuracy_metrics <- data.frame(
          TP = TP, FN = FN, FP = FP, TN = TN,
          sensitivity = sensitivity,
          specificity = specificity,
          PPV = PPV,
          NPV = NPV,
          proportion_area_covered = proportion_area_covered,
          jaccard_avg = jaccard_avg
        )
      }
      
      # Add state and level information
      accuracy_metrics$state <- state_name
      accuracy_metrics$level <- level
      accuracy_metrics$simulation <- sim
      
      # Store results
      sim_accuracy_metrics[[length(sim_accuracy_metrics) + 1]] <- accuracy_metrics
    } # End level loop
  } # End state loop
  
  # Combine accuracy metrics for this simulation
  if(length(sim_accuracy_metrics) > 0) {
    sim_df <- do.call(rbind, sim_accuracy_metrics)
    all_sim_results[[sim]] <- sim_df
  }
    
  }
} # End simulation loop

# Combine results from all simulations
if(length(all_sim_results) > 0) {
  final_results <- do.call(rbind, all_sim_results)
  
  # Calculate average metrics by state and level
  summary_results <- final_results %>%
    group_by(state, level) %>%
    summarise(
      avg_TP = mean(TP, na.rm = TRUE),
      avg_FN = mean(FN, na.rm = TRUE),
      avg_FP = mean(FP, na.rm = TRUE),
      avg_TN = mean(TN, na.rm = TRUE),
      avg_sensitivity = mean(sensitivity, na.rm = TRUE),
      avg_specificity = mean(specificity, na.rm = TRUE),
      avg_PPV = mean(PPV, na.rm = TRUE),
      avg_NPV = mean(NPV, na.rm = TRUE),
      avg_proportion_area_covered = mean(proportion_area_covered, na.rm = TRUE),
      avg_jaccard = mean(jaccard_avg, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Write final results to CSV
  write.csv(final_results, file.path(PATH_WD, "simulation_results_detailed.csv"), row.names = FALSE)
  write.csv(summary_results, file.path(PATH_WD, "simulation_results_summary.csv"), row.names = FALSE)
  
  cat("Simulations completed successfully!\n")
  cat("Detailed results saved to: simulation_results_detailed.csv\n")
  cat("Summary results saved to: simulation_results_summary.csv\n")
} else {
  cat("No results to save - all simulations failed.\n")
}
