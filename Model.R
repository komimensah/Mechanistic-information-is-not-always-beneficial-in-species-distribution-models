# ============================================================
# FINAL PUBLICATION PIPELINE
# Persistence-informed background selection for SDMs
# Species: Rhynchophorus phoenicis
# Domain: Global prediction at ~10 km
# Expert review: Africa-view blind maps
# Algorithms: GLM, GAM, RF, XGB, MaxEnt/maxnet, SVM
# Ensembles: mean ensemble + TSS-weighted ensemble
# ============================================================

rm(list = ls())
cat("\014")
set.seed(123)

# ============================================================
# 0. PACKAGES
# ============================================================

pkgs <- c(
  "terra", "sf", "dplyr", "mgcv", "ranger", "xgboost",
  "maxnet", "pROC", "ggplot2", "usdm", "corrplot",
  "e1071", "tidyterra", "viridis", "rnaturalearth",
  "rnaturalearthdata", "patchwork"
)

new_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if(length(new_pkgs) > 0) install.packages(new_pkgs)

library(terra)
library(sf)
library(dplyr)
library(mgcv)
library(ranger)
library(xgboost)
library(maxnet)
library(pROC)
library(ggplot2)
library(usdm)
library(corrplot)
library(e1071)
library(tidyterra)
library(viridis)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)

sf::sf_use_s2(FALSE)

# ============================================================
# 1. PATHS
# ============================================================

bio_dir  <- "/Users/kagboka/Downloads/wc2.1_30s_bio"
occ_path <- "/Users/kagboka/Desktop/Africa_palm_weevel/All Data.csv"

out_dir <- "/Users/kagboka/Desktop/Africa_palm_weevel/FINAL_SDM_background_comparison"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

species_name <- "Rhynchophorus_phoenicis"

# ============================================================
# 2. LOAD AND CLEAN OCCURRENCE DATA
# ============================================================

occ_raw <- read.csv(occ_path)
print(names(occ_raw))

lon_col <- names(occ_raw)[tolower(names(occ_raw)) %in% c(
  "lon", "longitude", "decimallongitude", "x"
)]

lat_col <- names(occ_raw)[tolower(names(occ_raw)) %in% c(
  "lat", "latitude", "decimallatitude", "y"
)]

if(length(lon_col) == 0 | length(lat_col) == 0){
  stop("Longitude/latitude columns not detected.")
}

occ_all <- occ_raw %>%
  dplyr::select(lon = all_of(lon_col[1]), lat = all_of(lat_col[1]), everything()) %>%
  distinct() %>%
  filter(!is.na(lon), !is.na(lat))

bad_coords <- occ_all %>%
  filter(lon < -180 | lon > 180 | lat < -90 | lat > 90)

write.csv(
  bad_coords,
  file.path(out_dir, "removed_invalid_coordinates.csv"),
  row.names = FALSE
)

occ_clean <- occ_all %>%
  filter(lon >= -180, lon <= 180, lat >= -90, lat <= 90)

world_sf <- rnaturalearth::ne_countries(returnclass = "sf")
africa_sf <- rnaturalearth::ne_countries(continent = "Africa", returnclass = "sf")

occ_sf <- st_as_sf(
  occ_clean,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
)

inside_africa <- st_intersects(occ_sf, africa_sf, sparse = FALSE)

occ_clean$in_africa <- rowSums(inside_africa) > 0

outside_africa <- occ_clean %>% filter(in_africa == FALSE)

write.csv(
  outside_africa,
  file.path(out_dir, "records_outside_africa_check_manually.csv"),
  row.names = FALSE
)

# Main analysis uses Africa records only because this is African palm weevil
occ <- occ_clean %>%
  filter(in_africa == TRUE) %>%
  dplyr::select(lon, lat) %>%
  distinct()

write.csv(
  occ,
  file.path(out_dir, "occurrences_used_for_modelling.csv"),
  row.names = FALSE
)

cat("Invalid coordinates removed:", nrow(bad_coords), "\n")
cat("Records outside Africa removed:", nrow(outside_africa), "\n")
cat("Final modelling records:", nrow(occ), "\n")

# ============================================================
# 3. LOAD WORLDCLIM AND AGGREGATE TO ~10 KM
# ============================================================

bio_files <- list.files(bio_dir, pattern = "\\.tif$", full.names = TRUE)

if(length(bio_files) == 0){
  stop("No .tif files found in bio_dir.")
}

bio_files_sorted <- bio_files[
  order(as.numeric(gsub(".*bio_?([0-9]+).*", "\\1", basename(bio_files))))
]

bio_1km <- rast(bio_files_sorted)
names(bio_1km) <- paste0("bio", 1:nlyr(bio_1km))

# 30 arc-sec ~1 km; aggregate by 10 = ~10 km
bio <- terra::aggregate(bio_1km, fact = 10, fun = mean, na.rm = TRUE)

# Global prediction domain
# Africa calibration domain
africa_vect <- terra::vect(africa_sf)

bio_crop <- terra::crop(bio, africa_vect)
bio_crop <- terra::mask(bio_crop, africa_vect)

# Keep global raster separately for final projection
bio_global <- bio
# ============================================================
# 4. VARIABLE SELECTION: CORRELATION DENDROGRAM + VIF
# ============================================================

var_dir <- file.path(out_dir, "variable_selection")
dir.create(var_dir, showWarnings = FALSE, recursive = TRUE)

bio_sample <- spatSample(
  bio_crop,
  size = 10000,
  method = "random",
  na.rm = TRUE,
  as.df = TRUE
)

bio_sample <- na.omit(bio_sample)

cor_mat <- cor(bio_sample, use = "complete.obs", method = "pearson")

write.csv(
  cor_mat,
  file.path(var_dir, "BIO_correlation_matrix.csv")
)

png(file.path(var_dir, "BIO_correlation_matrix.png"),
    width = 2200, height = 2000, res = 250)
corrplot(
  cor_mat,
  method = "color",
  type = "upper",
  tl.cex = 0.8,
  addCoef.col = "black",
  number.cex = 0.5
)
dev.off()

dist_mat <- as.dist(1 - abs(cor_mat))
hc <- hclust(dist_mat, method = "average")

png(file.path(var_dir, "BIO_correlation_dendrogram.png"),
    width = 1800, height = 1200, res = 250)
plot(hc, main = "BIOCLIM variable clustering", xlab = "", sub = "")
abline(h = 0.3, col = "red", lty = 2, lwd = 2)
dev.off()

clusters <- cutree(hc, h = 0.3)

cluster_df <- data.frame(
  variable = names(clusters),
  cluster = clusters
)

selected_corr <- c()

for(cl in unique(cluster_df$cluster)) {
  
  vars <- cluster_df$variable[cluster_df$cluster == cl]
  
  if(length(vars) == 1) {
    selected_corr <- c(selected_corr, vars)
  } else {
    sub_cor <- abs(cor_mat[vars, vars])
    mean_cor <- rowMeans(sub_cor, na.rm = TRUE)
    selected_corr <- c(selected_corr, names(which.min(mean_cor)))
  }
}

cat("Variables after dendrogram filtering:\n")
print(selected_corr)

bio_corr <- bio_crop[[selected_corr]]

vif_sample <- spatSample(
  bio_corr,
  size = 10000,
  method = "random",
  na.rm = TRUE,
  as.df = TRUE
)

vif_sample <- na.omit(vif_sample)

vif_res <- usdm::vifstep(vif_sample, th = 5)
selected_vars <- vif_res@results$Variables

cat("Final selected variables after VIF:\n")
print(selected_vars)

write.csv(
  data.frame(selected_variables = selected_vars),
  file.path(var_dir, "selected_BIO_variables.csv"),
  row.names = FALSE
)

bio_sel <- bio_crop[[selected_vars]]

bio1 <- bio_crop[["bio1"]]

# ============================================================
# 5. THERMAL FUNCTIONS AND PERSISTENCE METRIC
# ============================================================

egg.dev.rate <- function(T) {
  rep(0.244, length(T))
}

larva.dev.rate <- function(T) {
  b <- 0.0082
  Tmin <- 15.0316
  c <- 0.0364
  
  r <- ifelse(T >= Tmin, (b * (T - Tmin))^2 + c, 0)
  r[r < 0] <- 0
  r[r > 1] <- 1
  return(r)
}

pupa.dev.rate <- function(T) {
  a <- 1.609
  b <- 2.5613
  c <- 26.6194
  
  r <- a * (b / c) * (T / c)^(b - 1) * exp(-(T / c)^b)
  r[r < 0] <- 0
  r[r > 1] <- 1
  return(r)
}

egg.mort.rate <- function(T) {
  a <- 0.0328
  Tmin <- 49.6356
  Tmax <- 336.6055
  e <- 366.5798
  
  m <- ifelse(T <= Tmax, a * T * (T - Tmin) * sqrt(pmax(0, Tmax - T)) + e, NA)
  m[is.na(m) | m < 0] <- 0
  m[m > 1] <- 1
  return(m)
}

larva.mort.rate <- function(T) {
  a <- 1110.17
  b <- -128.7583
  c <- 4.8884
  d <- -0.0589
  
  m <- a + b*T + c*T^2 + d*T^3
  m[m < 0] <- 0
  m[m > 1] <- 1
  return(m)
}

pupa.mort.rate <- function(T) {
  a <- 4.5647
  b <- -0.2882
  c <- 0.0047
  
  m <- a + b*T + c*T^2
  m[m < 0] <- 0
  m[m > 1] <- 1
  return(m)
}

female.mort.rate <- function(T) {
  a <- 0.2494
  b <- -0.0207
  c <- 0.0005
  
  m <- a + b*T + c*T^2
  m[m < 0] <- 0
  m[m > 1] <- 1
  return(m)
}

fem.tot.ovip <- function(T) {
  a <- 154.9973
  c <- 24.723
  d <- 7.1077
  
  r <- a * exp(-((T - c) / d)^2)
  r[r < 0] <- 0
  return(r)
}

fem.ovip.rate <- function(T) {
  rep(6.12, length(T))
}

risk.index <- function(T) {
  
  sex.ratio <- 0.5
  
  Num <- egg.dev.rate(T) *
    larva.dev.rate(T) *
    pupa.dev.rate(T) *
    fem.ovip.rate(T) *
    fem.tot.ovip(T) *
    sex.ratio
  
  Den <- (egg.dev.rate(T) + egg.mort.rate(T)) *
    (larva.dev.rate(T) + larva.mort.rate(T)) *
    (pupa.dev.rate(T) + pupa.mort.rate(T)) *
    female.mort.rate(T)
  
  r <- Num / Den
  r[is.nan(r) | is.infinite(r) | r < 0] <- 0
  return(r)
}
# WorldClim BIO1 is usually °C 
bio1_temp <- bio1

risk_raster <- app(bio1_temp, fun = function(x) {
  r <- risk.index(x)
  r[!is.finite(r)] <- NA
  return(r)
})

names(risk_raster) <- "risk_index"

persistence <- app(bio1_temp, fun = function(x) {
  r <- risk.index(x)
  p <- ifelse(r > 1, 1 - (1 / r), 0)
  p[!is.finite(p)] <- NA
  return(p)
})

names(persistence) <- "persistence"

writeRaster(
  risk_raster,
  file.path(out_dir, "Rhynchophorus_phoenicis_raw_risk_index_global_10km.tif"),
  overwrite = TRUE
)

writeRaster(
  persistence,
  file.path(out_dir, "Rhynchophorus_phoenicis_branching_persistence_global_10km.tif"),
  overwrite = TRUE
)

predictors_base <- bio_sel
predictors_persist <- c(bio_sel, persistence)

# ============================================================
# 6. BACKGROUND SAMPLING FUNCTIONS
# ============================================================

sample_random_bg <- function(pred_stack, n = 10000) {
  
  bg <- spatSample(
    pred_stack[[1]],
    size = n,
    method = "random",
    na.rm = TRUE,
    as.points = TRUE
  )
  
  as.data.frame(crds(bg)) %>%
    setNames(c("lon", "lat"))
}

sample_kernel_bg <- function(pred_stack, occ_df, n = 10000) {
  
  template <- pred_stack[[1]]
  
  occ_points <- vect(
    occ_df,
    geom = c("lon", "lat"),
    crs = "EPSG:4326"
  )
  
  occ_r <- rasterize(
    occ_points,
    template,
    field = 1,
    fun = "sum",
    background = 0
  )
  
  k <- matrix(1, nrow = 21, ncol = 21)
  
  dens <- focal(
    occ_r,
    w = k,
    fun = sum,
    na.policy = "omit",
    fillvalue = 0
  )
  
  dens <- dens + 1e-6
  dens <- mask(dens, template)
  
  bg <- spatSample(
    dens,
    size = n,
    method = "weights",
    na.rm = TRUE,
    as.points = TRUE
  )
  
  as.data.frame(crds(bg)) %>%
    setNames(c("lon", "lat"))
}

sample_persistence_zero_bg <- function(pred_stack,
                                       persistence_raster,
                                       n = 10000) {
  
  template <- pred_stack[[1]]
  
  # Use the lowest 50% of persistence values
  threshold <- as.numeric(
    quantile(
      terra::values(persistence_raster),
      probs = 0.50,
      na.rm = TRUE
    )
  )
  
  message("Persistence threshold = ", round(threshold, 4))
  
  low_area <- persistence_raster
  
  # Keep only low-persistence cells
  low_area[low_area > threshold] <- NA
  low_area[low_area <= threshold] <- 1
  
  low_area <- terra::mask(low_area, template)
  
  available_cells <- terra::global(
    !is.na(low_area),
    "sum",
    na.rm = TRUE
  )[1,1]
  
  if(available_cells < n){
    warning("Fewer candidate cells than requested. Sampling all available cells.")
    n <- available_cells
  }
  
  bg <- terra::spatSample(
    low_area,
    size = n,
    method = "random",
    na.rm = TRUE,
    as.points = TRUE
  )
  
  bg_df <- as.data.frame(terra::crds(bg))
  names(bg_df) <- c("lon", "lat")
  
  return(bg_df)
}
# ============================================================
# 7. CREATE BACKGROUNDS
# ============================================================

n_bg <- 5000

bg_random <- sample_random_bg(
  predictors_base,
  n = n_bg
)

bg_kernel <- sample_kernel_bg(
  predictors_base,
  occ,
  n = n_bg
)

bg_persist_zero <- sample_persistence_zero_bg(
  predictors_base,
  persistence,
  n = n_bg
)
# ============================================================
# 8. MODELLING DATASETS
# ============================================================

make_dataset <- function(occ_df, bg_df, pred_stack) {
  
  pres_vals <- terra::extract(pred_stack, occ_df[, c("lon", "lat")])
  bg_vals   <- terra::extract(pred_stack, bg_df[, c("lon", "lat")])
  
  pres_dat <- cbind(pa = 1, pres_vals[, -1])
  bg_dat   <- cbind(pa = 0, bg_vals[, -1])
  
  dat <- rbind(pres_dat, bg_dat)
  dat <- na.omit(dat)
  dat$pa <- factor(dat$pa, levels = c(0, 1))
  
  return(dat)
}

# ============================================================
# 9. EVALUATION FUNCTION
# ============================================================

eval_model <- function(obs, pred) {
  
  obs_num <- as.numeric(as.character(obs))
  pred[!is.finite(pred)] <- NA
  
  keep <- complete.cases(obs_num, pred)
  obs_num <- obs_num[keep]
  pred <- pred[keep]
  
  auc_val <- as.numeric(pROC::auc(obs_num, pred, quiet = TRUE))
  
  thresholds <- seq(0.01, 0.99, by = 0.01)
  
  tss_values <- sapply(thresholds, function(th) {
    
    bin <- ifelse(pred >= th, 1, 0)
    
    TP <- sum(bin == 1 & obs_num == 1)
    TN <- sum(bin == 0 & obs_num == 0)
    FP <- sum(bin == 1 & obs_num == 0)
    FN <- sum(bin == 0 & obs_num == 1)
    
    sens <- ifelse((TP + FN) == 0, NA, TP / (TP + FN))
    spec <- ifelse((TN + FP) == 0, NA, TN / (TN + FP))
    
    sens + spec - 1
  })
  
  best_th <- thresholds[which.max(tss_values)]
  
  bin <- ifelse(pred >= best_th, 1, 0)
  
  TP <- sum(bin == 1 & obs_num == 1)
  TN <- sum(bin == 0 & obs_num == 0)
  FP <- sum(bin == 1 & obs_num == 0)
  FN <- sum(bin == 0 & obs_num == 1)
  
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  tss <- sensitivity + specificity - 1
  
  precision <- TP / (TP + FP)
  f1 <- 2 * precision * sensitivity / (precision + sensitivity)
  
  data.frame(
    AUC = auc_val,
    Threshold = best_th,
    Accuracy = accuracy,
    Sensitivity = sensitivity,
    Specificity = specificity,
    TSS = tss,
    F1 = f1
  )
}

# ============================================================
# 10. MODEL FUNCTIONS
# ============================================================

fit_glm <- function(train) {
  glm(pa ~ ., data = train, family = binomial)
}

# Replace fit_gam
fit_gam <- function(train) {
  
  vars <- names(train)[names(train) != "pa"]
  
  terms <- sapply(vars, function(v) {
    n_unique <- length(unique(train[[v]]))
    
    if(n_unique <= 5) {
      return(v)
    } else {
      k_use <- min(4, n_unique - 1)
      return(paste0("s(", v, ", k=", k_use, ")"))
    }
  })
  
  form <- as.formula(paste("pa ~", paste(terms, collapse = " + ")))
  
  mgcv::gam(form, data = train, family = binomial, method = "REML")
}

fit_rf <- function(train) {
  ranger::ranger(
    pa ~ .,
    data = train,
    probability = TRUE,
    num.trees = 500,
    importance = "impurity"
  )
}

# Replace fit_maxnet
fit_maxnet <- function(train) {
  
  x <- train[, names(train) != "pa", drop = FALSE]
  y <- as.numeric(as.character(train$pa))
  
  keep <- sapply(x, function(z) length(unique(z[is.finite(z)])) > 2)
  x <- x[, keep, drop = FALSE]
  
  maxnet::maxnet(
    p = y,
    data = x,
    f = maxnet::maxnet.formula(y, x, classes = "lqph")
  )
}

fit_svm <- function(train) {
  
  e1071::svm(
    pa ~ .,
    data = train,
    type = "C-classification",
    kernel = "radial",
    probability = TRUE,
    scale = TRUE
  )
}

predict_model <- function(model, newdata, algorithm) {
  
  newdata <- as.data.frame(newdata)
  
  if(algorithm == "GLM") {
    return(as.numeric(predict(model, newdata = newdata, type = "response")))
  }
  
  if(algorithm == "GAM") {
    return(as.numeric(predict(model, newdata = newdata, type = "response")))
  }
  
  if(algorithm == "RF") {
    return(as.numeric(predict(model, data = newdata)$predictions[, "1"]))
  }
  
  if(algorithm == "XGB") {
    return(as.numeric(predict(model, as.matrix(newdata))))
  }
  
  if(algorithm == "MAXNET") {
    needed <- names(model$samplemeans)
    needed <- needed[needed %in% names(newdata)]
    return(as.numeric(predict(model, newdata[, needed, drop = FALSE], type = "cloglog")))
  }
  
  if(algorithm == "SVM") {
    pr <- predict(model, newdata, probability = TRUE)
    prob <- attr(pr, "probabilities")
    
    if("1" %in% colnames(prob)) {
      return(as.numeric(prob[, "1"]))
    } else {
      return(as.numeric(prob[, ncol(prob)]))
    }
  }
}

fit_one_model <- function(dat, alg) {
  
  switch(
    alg,
    "GLM" = fit_glm(dat),
    "GAM" = fit_gam(dat),
    "RF" = fit_rf(dat),
    "XGB" = fit_xgb(dat),
    "MAXNET" = fit_maxnet(dat),
    "SVM" = fit_svm(dat)
  )
}

# ============================================================
# 11. CROSS-VALIDATION WITH ENSEMBLES
# ============================================================

create_folds <- function(dat, k = 5) {
  
  pres_id <- which(dat$pa == 1)
  abs_id  <- which(dat$pa == 0)
  
  pres_fold <- sample(rep(1:k, length.out = length(pres_id)))
  abs_fold  <- sample(rep(1:k, length.out = length(abs_id)))
  
  dat$fold <- NA
  dat$fold[pres_id] <- pres_fold
  dat$fold[abs_id] <- abs_fold
  
  dat
}

run_cv_with_ensembles <- function(dat, algorithms, k = 5) {
  
  dat <- create_folds(dat, k = k)
  out <- list()
  
  for(fold in 1:k) {
    
    train <- dat[dat$fold != fold, ]
    test  <- dat[dat$fold == fold, ]
    
    train$fold <- NULL
    test$fold <- NULL
    
    test_x <- test[, names(test) != "pa"]
    
    fold_preds <- list()
    fold_metrics <- list()
    
    for(alg in algorithms) {
      
      model <- fit_one_model(train, alg)
      pred <- predict_model(model, test_x, alg)
      
      fold_preds[[alg]] <- pred
      
      metrics <- eval_model(test$pa, pred)
      metrics$Fold <- fold
      metrics$Algorithm <- alg
      
      fold_metrics[[alg]] <- metrics
    }
    
    pred_mat <- do.call(cbind, fold_preds)
    
    ens_mean <- rowMeans(pred_mat, na.rm = TRUE)
    
    m_mean <- eval_model(test$pa, ens_mean)
    m_mean$Fold <- fold
    m_mean$Algorithm <- "ENSEMBLE_MEAN"
    
    indiv_tss <- sapply(fold_metrics, function(x) x$TSS[1])
    indiv_tss[indiv_tss < 0 | !is.finite(indiv_tss)] <- 0
    
    if(sum(indiv_tss, na.rm = TRUE) == 0) {
      weights <- rep(1 / length(indiv_tss), length(indiv_tss))
    } else {
      weights <- indiv_tss / sum(indiv_tss, na.rm = TRUE)
    }
    
    ens_weighted <- as.numeric(pred_mat %*% weights)
    
    m_weighted <- eval_model(test$pa, ens_weighted)
    m_weighted$Fold <- fold
    m_weighted$Algorithm <- "ENSEMBLE_WEIGHTED_TSS"
    
    out[[paste0("fold_", fold)]] <- bind_rows(
      bind_rows(fold_metrics),
      m_mean,
      m_weighted
    )
  }
  
  bind_rows(out)
}

# ============================================================
# 12. TREATMENTS
# ============================================================

treatments <- list(
  Random_background = list(
    bg = bg_random,
    predictors = predictors_base
  ),
  Kernel_background = list(
    bg = bg_kernel,
    predictors = predictors_base
  ),
  Persistence_zero_background = list(
    bg = bg_persist_zero,
    predictors = predictors_base
  ),
  Persistence_predictor = list(
    bg = bg_random,
    predictors = predictors_persist
  ),
  
  # ADD THIS NEW TREATMENT HERE
  Kernel_background_plus_predictor = list(
    bg = bg_kernel,
    predictors = predictors_persist
  ),
  
  Persistence_zero_background_plus_predictor = list(
    bg = bg_persist_zero,
    predictors = predictors_persist
  )
)

algorithms <- c("GLM", "GAM", "RF", "MAXNET", "SVM")

# ============================================================
# 13. RUN CROSS-VALIDATION -- SAFE VERSION
# ============================================================

all_results <- list()

for(tr in names(treatments)) {
  
  message("Running treatment: ", tr)
  
  dat <- make_dataset(
    occ,
    treatments[[tr]]$bg,
    treatments[[tr]]$predictors
  )
  
  message("  Dataset rows: ", nrow(dat))
  message("  Presence: ", sum(dat$pa == 1))
  message("  Background: ", sum(dat$pa == 0))
  
  cv_metrics <- tryCatch({
    run_cv_with_ensembles(dat, algorithms, k = 5)
  }, error = function(e) {
    message("ERROR in ", tr, ": ", e$message)
    return(NULL)
  })
  
  if(!is.null(cv_metrics) && nrow(cv_metrics) > 0) {
    
    cv_metrics$Treatment <- tr
    cv_metrics$N_presence <- sum(dat$pa == 1)
    cv_metrics$N_background <- sum(dat$pa == 0)
    
    all_results[[tr]] <- cv_metrics
  }
}

if(length(all_results) == 0) {
  stop("No cross-validation results were produced. Check the error messages above.")
}

results_cv <- bind_rows(all_results)

print(names(results_cv))
print(head(results_cv))

write.csv(
  results_cv,
  file.path(out_dir, "cross_validation_results_with_ensembles.csv"),
  row.names = FALSE
)

summary_results <- results_cv %>%
  group_by(Treatment, Algorithm) %>%
  summarise(
    AUC_mean = mean(AUC, na.rm = TRUE),
    AUC_sd = sd(AUC, na.rm = TRUE),
    TSS_mean = mean(TSS, na.rm = TRUE),
    TSS_sd = sd(TSS, na.rm = TRUE),
    F1_mean = mean(F1, na.rm = TRUE),
    F1_sd = sd(F1, na.rm = TRUE),
    Accuracy_mean = mean(Accuracy, na.rm = TRUE),
    Sensitivity_mean = mean(Sensitivity, na.rm = TRUE),
    Specificity_mean = mean(Specificity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(TSS_mean))

write.csv(
  summary_results,
  file.path(out_dir, "cross_validation_summary_with_ensembles.csv"),
  row.names = FALSE
)

print(summary_results)
# ============================================================
# 13B. SENSITIVITY ANALYSIS:
# Persistence-background threshold gradient
# 5% to 50% lower persistence values
# ============================================================

sensitivity_dir <- file.path(out_dir, "sensitivity_persistence_threshold_gradient")
dir.create(sensitivity_dir, showWarnings = FALSE, recursive = TRUE)

sample_persistence_threshold_bg <- function(pred_stack,
                                            persistence_raster,
                                            n = 5000,
                                            threshold_prob = 0.50) {
  
  template <- pred_stack[[1]]
  
  threshold <- as.numeric(
    quantile(
      terra::values(persistence_raster),
      probs = threshold_prob,
      na.rm = TRUE
    )
  )
  
  low_area <- persistence_raster
  low_area[low_area > threshold] <- NA
  low_area[low_area <= threshold] <- 1
  low_area <- terra::mask(low_area, template)
  
  available_cells <- terra::global(!is.na(low_area), "sum", na.rm = TRUE)[1, 1]
  
  if(available_cells < n) {
    warning("Fewer candidate cells than requested. Sampling all available cells.")
    n <- available_cells
  }
  
  bg <- terra::spatSample(
    low_area,
    size = n,
    method = "random",
    na.rm = TRUE,
    as.points = TRUE
  )
  
  bg_df <- as.data.frame(terra::crds(bg))
  names(bg_df) <- c("lon", "lat")
  
  return(bg_df)
}

# Threshold gradient: from very strict to more relaxed
threshold_probs <- seq(0.05, 0.50, by = 0.05)

sensitivity_results <- list()

for(th_prob in threshold_probs) {
  
  message("Running persistence-threshold sensitivity: ", th_prob)
  
  bg_persist_sens <- sample_persistence_threshold_bg(
    pred_stack = predictors_base,
    persistence_raster = persistence,
    n = n_bg,
    threshold_prob = th_prob
  )
  
  treatments_sens <- list(
    Persistence_background = list(
      bg = bg_persist_sens,
      predictors = predictors_base
    ),
    Persistence_background_plus_predictor = list(
      bg = bg_persist_sens,
      predictors = predictors_persist
    )
  )
  
  for(tr in names(treatments_sens)) {
    
    message("  Treatment: ", tr)
    
    dat <- make_dataset(
      occ,
      treatments_sens[[tr]]$bg,
      treatments_sens[[tr]]$predictors
    )
    
    cv_metrics <- tryCatch({
      run_cv_with_ensembles(dat, algorithms, k = 5)
    }, error = function(e) {
      message("ERROR in threshold ", th_prob, " treatment ", tr, ": ", e$message)
      return(NULL)
    })
    
    if(!is.null(cv_metrics)) {
      cv_metrics$Treatment <- tr
      cv_metrics$Threshold_prob <- th_prob
      cv_metrics$Threshold_percent <- th_prob * 100
      sensitivity_results[[paste0(tr, "_", th_prob)]] <- cv_metrics
    }
  }
}

sensitivity_cv <- dplyr::bind_rows(sensitivity_results)

write.csv(
  sensitivity_cv,
  file.path(sensitivity_dir, "threshold_gradient_cross_validation_raw.csv"),
  row.names = FALSE
)

# ============================================================
# Summarise with 95% confidence intervals
# ============================================================

sensitivity_summary <- sensitivity_cv %>%
  group_by(Treatment, Threshold_prob, Threshold_percent) %>%
  summarise(
    AUC_mean = mean(AUC, na.rm = TRUE),
    AUC_sd = sd(AUC, na.rm = TRUE),
    AUC_n = sum(!is.na(AUC)),
    AUC_se = AUC_sd / sqrt(AUC_n),
    AUC_lwr = AUC_mean - 1.96 * AUC_se,
    AUC_upr = AUC_mean + 1.96 * AUC_se,
    
    TSS_mean = mean(TSS, na.rm = TRUE),
    TSS_sd = sd(TSS, na.rm = TRUE),
    TSS_n = sum(!is.na(TSS)),
    TSS_se = TSS_sd / sqrt(TSS_n),
    TSS_lwr = TSS_mean - 1.96 * TSS_se,
    TSS_upr = TSS_mean + 1.96 * TSS_se,
    
    .groups = "drop"
  )

write.csv(
  sensitivity_summary,
  file.path(sensitivity_dir, "threshold_gradient_cross_validation_summary_CI.csv"),
  row.names = FALSE
)

print(sensitivity_summary)

# ============================================================
# Plot TSS sensitivity
# ============================================================

library(ggplot2)
library(patchwork)

treatment_labels_sens <- c(
  "Persistence_background" = "Persistence background",
  "Persistence_background_plus_predictor" = "Persistence background + predictor"
)

plot_sens <- sensitivity_summary %>%
  mutate(
    Treatment_label = factor(
      treatment_labels_sens[Treatment],
      levels = c(
        "Persistence background",
        "Persistence background + predictor"
      )
    )
  )

p_tss_sens <- ggplot(
  plot_sens,
  aes(
    x = Threshold_percent,
    y = TSS_mean,
    ymin = TSS_lwr,
    ymax = TSS_upr,
    group = Treatment_label,
    linetype = Treatment_label
  )
) +
  geom_ribbon(alpha = 0.18, aes(fill = Treatment_label), color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_bw(base_size = 12) +
  labs(
    x = "Persistence threshold used for background sampling (%)",
    y = "Mean TSS",
    linetype = NULL,
    fill = NULL,
    title = "Sensitivity of TSS to persistence-background threshold"
  ) +
  scale_x_continuous(breaks = threshold_probs * 100) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top"
  )

# ============================================================
# Plot AUC sensitivity
# ============================================================

p_auc_sens <- ggplot(
  plot_sens,
  aes(
    x = Threshold_percent,
    y = AUC_mean,
    ymin = AUC_lwr,
    ymax = AUC_upr,
    group = Treatment_label,
    linetype = Treatment_label
  )
) +
  geom_ribbon(alpha = 0.18, aes(fill = Treatment_label), color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_bw(base_size = 12) +
  labs(
    x = "Persistence threshold used for background sampling (%)",
    y = "Mean AUC",
    linetype = NULL,
    fill = NULL,
    title = "Sensitivity of AUC to persistence-background threshold"
  ) +
  scale_x_continuous(breaks = threshold_probs * 100) +
  coord_cartesian(ylim = c(0.5, 1)) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top"
  )

fig_sensitivity <- p_tss_sens / p_auc_sens

ggsave(
  file.path(sensitivity_dir, "Supplementary_Figure_persistence_threshold_sensitivity_CI.png"),
  fig_sensitivity,
  width = 10,
  height = 9,
  dpi = 600
)

ggsave(
  file.path(sensitivity_dir, "Supplementary_Figure_persistence_threshold_sensitivity_CI.pdf"),
  fig_sensitivity,
  width = 10,
  height = 9
)

print(fig_sensitivity)

# ============================================================
# Manuscript-quality sensitivity figure with 95% CI ribbons
# ============================================================

library(ggplot2)
library(patchwork)
library(dplyr)

treatment_labels_sens <- c(
  "Persistence_background" = "Persistence background",
  "Persistence_background_plus_predictor" = "Persistence background + predictor"
)

plot_sens <- sensitivity_summary %>%
  mutate(
    Treatment_label = factor(
      treatment_labels_sens[Treatment],
      levels = c(
        "Persistence background",
        "Persistence background + predictor"
      )
    )
  )

# Optional: keep CI inside valid metric range
plot_sens <- plot_sens %>%
  mutate(
    TSS_lwr = pmax(TSS_lwr, 0),
    TSS_upr = pmin(TSS_upr, 1),
    AUC_lwr = pmax(AUC_lwr, 0.5),
    AUC_upr = pmin(AUC_upr, 1)
  )

p_tss_sens <- ggplot(
  plot_sens,
  aes(
    x = Threshold_percent,
    y = TSS_mean,
    group = Treatment_label
  )
) +
  geom_ribbon(
    aes(ymin = TSS_lwr, ymax = TSS_upr, fill = Treatment_label),
    alpha = 0.22,
    colour = NA
  ) +
  geom_line(
    aes(linetype = Treatment_label),
    linewidth = 1.1
  ) +
  geom_point(size = 2.4) +
  scale_x_continuous(
    breaks = seq(5, 50, 5),
    limits = c(5, 50)
  ) +
  coord_cartesian(ylim = c(0.80, 1.00)) +
  theme_bw(base_size = 12) +
  labs(
    x = NULL,
    y = "Mean TSS",
    fill = NULL,
    linetype = NULL,
    title = "A. TSS sensitivity"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25),
    plot.title = element_text(face = "bold", hjust = 0),
    legend.position = "top",
    legend.box = "horizontal"
  )

p_auc_sens <- ggplot(
  plot_sens,
  aes(
    x = Threshold_percent,
    y = AUC_mean,
    group = Treatment_label
  )
) +
  geom_ribbon(
    aes(ymin = AUC_lwr, ymax = AUC_upr, fill = Treatment_label),
    alpha = 0.22,
    colour = NA
  ) +
  geom_line(
    aes(linetype = Treatment_label),
    linewidth = 1.1
  ) +
  geom_point(size = 2.4) +
  scale_x_continuous(
    breaks = seq(5, 50, 5),
    limits = c(5, 50)
  ) +
  coord_cartesian(ylim = c(0.90, 1.00)) +
  theme_bw(base_size = 12) +
  labs(
    x = "Persistence lower-tail threshold used for background sampling (%)",
    y = "Mean AUC",
    fill = NULL,
    linetype = NULL,
    title = "B. AUC sensitivity"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25),
    plot.title = element_text(face = "bold", hjust = 0),
    legend.position = "none"
  )

fig_sensitivity <- p_tss_sens / p_auc_sens +
  plot_annotation(
    title = "Sensitivity of persistence-informed background calibration to threshold choice",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5)
    )
  )

ggsave(
  file.path(sensitivity_dir, "Supplementary_Figure_threshold_sensitivity_curves_CI.png"),
  fig_sensitivity,
  width = 10,
  height = 8,
  dpi = 600
)

ggsave(
  file.path(sensitivity_dir, "Supplementary_Figure_threshold_sensitivity_curves_CI.pdf"),
  fig_sensitivity,
  width = 10,
  height = 8
)

print(fig_sensitivity)
# ============================================================
# Better sensitivity figure:
# Boxplots + jittered model/fold values + mean trend
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)

sensitivity_fig_dir <- file.path(sensitivity_dir, "better_sensitivity_figures")
dir.create(sensitivity_fig_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# Prepare long-format data
# ------------------------------------------------------------

sens_long <- sensitivity_cv %>%
  dplyr::select(
    Treatment,
    Algorithm,
    Fold,
    Threshold_percent,
    AUC,
    TSS
  ) %>%
  tidyr::pivot_longer(
    cols = c(AUC, TSS),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Treatment_label = dplyr::recode(
      Treatment,
      "Persistence_background" = "Persistence background",
      "Persistence_background_plus_predictor" = "Persistence background + predictor"
    ),
    Metric = factor(Metric, levels = c("TSS", "AUC")),
    Threshold_percent = factor(
      Threshold_percent,
      levels = sort(unique(Threshold_percent))
    )
  )

# ------------------------------------------------------------
# Mean and 95% CI by threshold/treatment/metric
# ------------------------------------------------------------

sens_mean <- sens_long %>%
  group_by(Treatment_label, Metric, Threshold_percent) %>%
  summarise(
    mean_value = mean(Value, na.rm = TRUE),
    sd_value = sd(Value, na.rm = TRUE),
    n = sum(!is.na(Value)),
    se = sd_value / sqrt(n),
    lwr = mean_value - 1.96 * se,
    upr = mean_value + 1.96 * se,
    .groups = "drop"
  ) %>%
  mutate(
    Threshold_numeric = as.numeric(as.character(Threshold_percent))
  )

# ------------------------------------------------------------
# Manuscript sensitivity figure
# ------------------------------------------------------------

p_sens_box <- ggplot(
  sens_long,
  aes(
    x = Threshold_percent,
    y = Value,
    fill = Treatment_label
  )
) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.6,
    alpha = 0.55,
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(color = Treatment_label),
    position = position_jitterdodge(
      jitter.width = 0.12,
      dodge.width = 0.75
    ),
    size = 1.2,
    alpha = 0.35
  ) +
  geom_line(
    data = sens_mean,
    aes(
      x = Threshold_numeric,
      y = mean_value,
      group = Treatment_label,
      color = Treatment_label
    ),
    inherit.aes = FALSE,
    linewidth = 1.1
  ) +
  geom_point(
    data = sens_mean,
    aes(
      x = Threshold_numeric,
      y = mean_value,
      color = Treatment_label
    ),
    inherit.aes = FALSE,
    size = 2.4
  ) +
  geom_errorbar(
    data = sens_mean,
    aes(
      x = Threshold_numeric,
      ymin = lwr,
      ymax = upr,
      color = Treatment_label
    ),
    inherit.aes = FALSE,
    width = 0.8,
    linewidth = 0.5
  ) +
  facet_wrap(
    ~ Metric,
    ncol = 1,
    scales = "free_y"
  ) +
  theme_bw(base_size = 12) +
  labs(
    x = "Lower-tail persistence threshold used for background sampling (%)",
    y = "Cross-validation performance",
    fill = NULL,
    color = NULL,
    title = "Sensitivity of persistence-informed background calibration to threshold choice"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "top",
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

ggsave(
  file.path(sensitivity_fig_dir, "Supplementary_Figure_threshold_sensitivity_boxplot_CI.png"),
  p_sens_box,
  width = 10,
  height = 8,
  dpi = 600
)

ggsave(
  file.path(sensitivity_fig_dir, "Supplementary_Figure_threshold_sensitivity_boxplot_CI.pdf"),
  p_sens_box,
  width = 10,
  height = 8
)

print(p_sens_box)


# ============================================================
# 14. RASTER PREDICTION -- GLOBAL PROJECTION
# Train on Africa, project globally
# ============================================================

# ------------------------------------------------------------
# 14.1 Build global predictor stacks
# ------------------------------------------------------------

bio_sel_global <- bio_global[[selected_vars]]
bio1_global <- bio_global[["bio1"]]

risk_raster_global <- app(bio1_global, fun = function(x) {
  r <- risk.index(x)
  r[!is.finite(r)] <- NA
  return(r)
})

names(risk_raster_global) <- "risk_index"

persistence_global <- app(bio1_global, fun = function(x) {
  r <- risk.index(x)
  p <- ifelse(r > 1, 1 - (1 / r), 0)
  p[!is.finite(p)] <- NA
  return(p)
})

names(persistence_global) <- "persistence"

predictors_base_global <- bio_sel_global
predictors_persist_global <- c(bio_sel_global, persistence_global)

# ------------------------------------------------------------
# 14.2 Prediction function
# ------------------------------------------------------------

predict_raster_model <- function(model, pred_stack, algorithm, filename) {
  
  terra::predict(
    pred_stack,
    model = model,
    fun = function(model, x) {
      
      x <- as.data.frame(x)
      
      if(ncol(x) == length(names(pred_stack))) {
        names(x) <- names(pred_stack)
      }
      
      out <- rep(NA_real_, nrow(x))
      keep <- stats::complete.cases(x)
      
      if(sum(keep) > 0) {
        pred <- tryCatch({
          predict_model(model, x[keep, , drop = FALSE], algorithm)
        }, error = function(e) {
          rep(NA_real_, sum(keep))
        })
        
        out[keep] <- as.numeric(pred)
      }
      
      return(out)
    },
    filename = filename,
    overwrite = TRUE,
    na.rm = FALSE
  )
}

map_dir <- file.path(out_dir, "prediction_maps_global_10km")
dir.create(map_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 14.3 Match Africa training stacks with global projection stacks
# ------------------------------------------------------------

projection_stacks <- list(
  Random_background = predictors_base_global,
  Kernel_background = predictors_base_global,
  Persistence_zero_background = predictors_base_global,
  Persistence_predictor = predictors_persist_global,
  Kernel_background_plus_predictor = predictors_persist_global,
  Persistence_zero_background_plus_predictor = predictors_persist_global
)

# ------------------------------------------------------------
# 14.4 Train on Africa and project globally
# ------------------------------------------------------------

for(tr in names(treatments)) {
  
  message("Mapping treatment globally: ", tr)
  
  # Africa stack for model training
  train_stack <- treatments[[tr]]$predictors
  
  # Global stack for projection
  pred_stack_global <- projection_stacks[[tr]]
  
  bg_df <- treatments[[tr]]$bg
  
  dat <- make_dataset(
    occ,
    bg_df,
    train_stack
  )
  
  model_map_files <- c()
  
  for(alg in algorithms) {
    
    message("  Mapping algorithm: ", alg)
    
    model <- tryCatch({
      fit_one_model(dat, alg)
    }, error = function(e) {
      message("Model failed: ", tr, " - ", alg, ": ", e$message)
      return(NULL)
    })
    
    if(!is.null(model)) {
      
      tif_file <- file.path(
        map_dir,
        paste0(species_name, "_", tr, "_", alg, ".tif")
      )
      
      pred_map <- tryCatch({
        predict_raster_model(
          model,
          pred_stack_global,
          alg,
          tif_file
        )
      }, error = function(e) {
        message("Raster prediction failed: ", tr, " - ", alg, ": ", e$message)
        return(NULL)
      })
      
      if(!is.null(pred_map) && file.exists(tif_file)) {
        model_map_files[alg] <- tif_file
      }
    }
  }
  
  # --------------------------
  # Global ensemble maps
  # --------------------------
  if(length(model_map_files) >= 2) {
    
    valid_algs <- names(model_map_files)
    
    stack_maps <- terra::rast(as.character(model_map_files))
    names(stack_maps) <- valid_algs
    
    ens_mean <- terra::app(
      stack_maps,
      fun = function(x) mean(x, na.rm = TRUE)
    )
    
    names(ens_mean) <- "ENSEMBLE_MEAN"
    
    writeRaster(
      ens_mean,
      file.path(map_dir, paste0(species_name, "_", tr, "_ENSEMBLE_MEAN.tif")),
      overwrite = TRUE
    )
    
    weights_df <- summary_results %>%
      dplyr::filter(Treatment == tr, Algorithm %in% valid_algs) %>%
      dplyr::select(Algorithm, TSS_mean)
    
    weights <- weights_df$TSS_mean
    names(weights) <- weights_df$Algorithm
    weights[weights < 0 | !is.finite(weights)] <- 0
    
    common_algs <- intersect(valid_algs, names(weights))
    
    if(length(common_algs) >= 2) {
      
      weights <- weights[common_algs]
      
      if(sum(weights) == 0) {
        weights <- rep(1 / length(weights), length(weights))
      } else {
        weights <- weights / sum(weights)
      }
      
      maps_weighted <- terra::rast(as.character(model_map_files[common_algs]))
      names(maps_weighted) <- common_algs
      
      ens_weighted <- terra::app(
        maps_weighted,
        fun = function(x) sum(x * weights, na.rm = TRUE)
      )
      
      names(ens_weighted) <- "ENSEMBLE_WEIGHTED_TSS"
      
      writeRaster(
        ens_weighted,
        file.path(map_dir, paste0(species_name, "_", tr, "_ENSEMBLE_WEIGHTED_TSS.tif")),
        overwrite = TRUE
      )
    }
  }
}

message("Global raster prediction completed. Maps saved in: ", map_dir)
# ============================================================
# 15. HIGH-QUALITY PUBLICATION AND EXPERT MAPS
# ============================================================

world <- rnaturalearth::ne_countries(returnclass = "sf")
africa <- rnaturalearth::ne_countries(continent = "Africa", returnclass = "sf")

plot_one_model_map <- function(tif_file, code = NULL, region = "Africa") {
  
  r <- rast(tif_file)
  nm <- tools::file_path_sans_ext(basename(tif_file))
  title_use <- ifelse(is.null(code), nm, code)
  
  if(region == "Africa") {
    boundary <- africa
    xlims <- c(-20, 55)
    ylims <- c(-35, 38)
  } else {
    boundary <- world
    xlims <- c(-180, 180)
    ylims <- c(-60, 80)
  }
  
  ggplot() +
    geom_sf(data = boundary, fill = "grey95", color = "grey35", linewidth = 0.2) +
    tidyterra::geom_spatraster(data = r) +
    geom_point(
      data = occ,
      aes(x = lon, y = lat),
      size = 0.35,
      alpha = 0.55
    ) +
    scale_fill_viridis_c(
      name = "Suitability",
      option = "C",
      limits = c(0, 1),
      na.value = "transparent"
    ) +
    coord_sf(xlim = xlims, ylim = ylims, expand = FALSE) +
    labs(title = title_use, x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      axis.text = element_text(color = "black", size = 7),
      plot.title = element_text(face = "bold", size = 9),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 8),
      legend.text = element_text(size = 7),
      legend.key.height = unit(0.9, "cm")
    )
}

tif_files <- list.files(map_dir, pattern = "\\.tif$", full.names = TRUE)

pub_global_dir <- file.path(out_dir, "publication_maps_global")
pub_africa_dir <- file.path(out_dir, "publication_maps_africa")
expert_dir <- file.path(out_dir, "expert_blind_maps_Africa")

dir.create(pub_global_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(pub_africa_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(expert_dir, showWarnings = FALSE, recursive = TRUE)

for(f in tif_files) {
  
  nm <- tools::file_path_sans_ext(basename(f))
  
  p_global <- plot_one_model_map(f, region = "Global")
  
  ggsave(
    file.path(pub_global_dir, paste0(nm, "_global.png")),
    p_global,
    width = 11,
    height = 6,
    dpi = 600
  )
  
  p_africa <- plot_one_model_map(f, region = "Africa")
  
  ggsave(
    file.path(pub_africa_dir, paste0(nm, "_africa.png")),
    p_africa,
    width = 7,
    height = 6,
    dpi = 600
  )
}

set.seed(999)
tif_files_random <- sample(tif_files)

blind_key <- data.frame(
  Code = paste0("Map_", sprintf("%02d", seq_along(tif_files_random))),
  Original_file = basename(tif_files_random)
)

for(i in seq_along(tif_files_random)) {
  
  code <- blind_key$Code[i]
  p <- plot_one_model_map(tif_files_random[i], code = code, region = "Africa")
  
  ggsave(
    file.path(expert_dir, paste0(code, ".png")),
    p,
    width = 7,
    height = 6,
    dpi = 600
  )
}

write.csv(
  blind_key,
  file.path(out_dir, "blind_map_key_DO_NOT_SEND_TO_EXPERTS.csv"),
  row.names = FALSE
)

expert_template <- data.frame(
  Code = blind_key$Code,
  Ecological_realism_1_to_5 = NA,
  Known_distribution_agreement_1_to_5 = NA,
  Overprediction_risk_1_to_5 = NA,
  Underprediction_risk_1_to_5 = NA,
  Management_usefulness_1_to_5 = NA,
  Overall_rank = NA,
  Comments = NA
)

write.csv(
  expert_template,
  file.path(out_dir, "expert_scoring_template.csv"),
  row.names = FALSE
)

# ============================================================
# 16. PERFORMANCE FIGURES
# ============================================================

p1 <- ggplot(
  summary_results,
  aes(x = Treatment, y = TSS_mean, fill = Algorithm)
) +
  geom_col(position = "dodge") +
  theme_bw() +
  coord_flip() +
  labs(
    x = "Treatment",
    y = "Mean TSS",
    title = "Cross-validation performance including ensembles"
  )

ggsave(
  file.path(out_dir, "Figure_TSS_summary_with_ensembles.png"),
  p1,
  width = 11,
  height = 8,
  dpi = 600
)

p2 <- ggplot(
  summary_results,
  aes(x = Treatment, y = AUC_mean, fill = Algorithm)
) +
  geom_col(position = "dodge") +
  theme_bw() +
  coord_flip() +
  labs(
    x = "Treatment",
    y = "Mean AUC",
    title = "Cross-validation AUC including ensembles"
  )

ggsave(
  file.path(out_dir, "Figure_AUC_summary_with_ensembles.png"),
  p2,
  width = 11,
  height = 8,
  dpi = 600
)

message("DONE. Final publication pipeline completed. Outputs saved in: ", out_dir)
# ============================================================
# HIGH-QUALITY SIDE-BY-SIDE TREATMENT MAPS
# Same style as previous grasshopper figure
# Creates Africa and Global figures for each algorithm
# ============================================================

library(terra)
library(ggplot2)
library(dplyr)
library(patchwork)
library(viridis)

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

map_dir <- file.path(out_dir, "prediction_maps_global_10km")

fig_dir <- file.path(out_dir, "publication_side_by_side_maps")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# Treatments and labels
# ------------------------------------------------------------

treatments_order <- c(
  "Random_background",
  "Kernel_background",
  "Persistence_zero_background",
  "Persistence_predictor",
  "Kernel_background_plus_predictor",
  "Persistence_zero_background_plus_predictor"
)

treatment_labels <- c(
  "Random_background" = "Random background",
  "Kernel_background" = "Kernel-density background",
  "Persistence_zero_background" = "Persistence-zero background",
  "Persistence_predictor" = "Persistence as predictor",
  "Kernel_background_plus_predictor" = "Kernel background + persistence predictor",
  "Persistence_zero_background_plus_predictor" = "Persistence background + predictor"
)

algorithm_labels <- c(
  "GLM" = "GLM",
  "GAM" = "GAM",
  "RF" = "Random Forest",
  "MAXNET" = "MaxEnt",
  "SVM" = "SVM",
  "ENSEMBLE_MEAN" = "Mean ensemble",
  "ENSEMBLE_WEIGHTED_TSS" = "TSS-weighted ensemble"
)

# ------------------------------------------------------------
# Function to read treatment rasters for one algorithm
# ------------------------------------------------------------

get_algorithm_stack <- function(algorithm_name) {
  
  files <- file.path(
    map_dir,
    paste0(species_name, "_", treatments_order, "_", algorithm_name, ".tif")
  )
  
  keep <- file.exists(files)
  
  if(sum(keep) == 0) {
    message("No files found for ", algorithm_name)
    return(NULL)
  }
  
  files <- files[keep]
  kept_treatments <- treatments_order[keep]
  
  r <- terra::rast(files)
  names(r) <- treatment_labels[kept_treatments]
  
  return(r)
}

# ------------------------------------------------------------
# Crop functions
# ------------------------------------------------------------

crop_africa_extent <- function(r) {
  terra::crop(r, terra::ext(-20, 55, -35, 38))
}

crop_global_extent <- function(r) {
  terra::crop(r, terra::ext(-180, 180, -60, 80))
}

# ------------------------------------------------------------
# Reduce resolution for plotting only
# ------------------------------------------------------------

reduce_for_plotting <- function(r, fact = 2) {
  
  if(fact > 1) {
    r <- terra::aggregate(
      r,
      fact = fact,
      fun = mean,
      na.rm = TRUE
    )
  }
  
  return(r)
}

# ------------------------------------------------------------
# Plot one raster layer
# ------------------------------------------------------------

plot_raster_clean <- function(r, title_text, show_legend = TRUE) {
  
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df) <- c("Longitude", "Latitude", "Suitability")
  
  p <- ggplot(df, aes(x = Longitude, y = Latitude, fill = Suitability)) +
    geom_raster() +
    coord_equal(expand = FALSE) +
    scale_fill_viridis(
      option = "magma",
      limits = c(0, 1),
      name = NULL,
      guide = guide_colorbar(
        barheight = unit(1.8, "cm"),
        barwidth  = unit(0.25, "cm")
      )
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 6),
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      legend.position = ifelse(show_legend, "right", "none"),
      legend.text = element_text(size = 6),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin = margin(2, 2, 2, 2)
    ) +
    labs(
      title = title_text,
      x = "Longitude",
      y = "Latitude"
    )
  
  return(p)
}

# ------------------------------------------------------------
# Build figure for one algorithm
# ------------------------------------------------------------

make_side_by_side_figure <- function(
    algorithm_name,
    region = "Africa",
    plot_fact = 2
) {
  
  r_stack <- get_algorithm_stack(algorithm_name)
  
  if(is.null(r_stack)) return(NULL)
  
  if(region == "Africa") {
    r_stack <- crop_africa_extent(r_stack)
  }
  
  if(region == "Global") {
    r_stack <- crop_global_extent(r_stack)
  }
  
  r_stack <- reduce_for_plotting(r_stack, fact = plot_fact)
  
  n_layers <- terra::nlyr(r_stack)
  
  plots <- list()
  
  for(i in seq_len(n_layers)) {
    
    panel_label <- paste0(letters[i], ") ")
    title_text <- paste0(panel_label, names(r_stack)[i])
    
    plots[[i]] <- plot_raster_clean(
      r_stack[[i]],
      title_text = title_text,
      show_legend = TRUE
    )
  }
  
  # layout: 5 panels = 3 top, 2 bottom
  if(n_layers == 5) {
    fig <- (plots[[1]] | plots[[2]] | plots[[3]]) /
      (plots[[4]] | plots[[5]] | patchwork::plot_spacer())
  } else if(n_layers == 4) {
    fig <- (plots[[1]] | plots[[2]]) /
      (plots[[3]] | plots[[4]])
  } else if(n_layers == 3) {
    fig <- plots[[1]] | plots[[2]] | plots[[3]]
  } else {
    fig <- patchwork::wrap_plots(plots)
  }
  
  fig <- fig +
    plot_annotation(
      title = paste0(
        algorithm_labels[[algorithm_name]],
        " suitability maps across background-selection treatments — ",
        region
      ),
      theme = theme(
        plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
      )
    )
  
  return(fig)
}

# ------------------------------------------------------------
# Create Africa and Global figures for all algorithms
# ------------------------------------------------------------

algorithms_to_plot <- c(
  "GLM",
  "GAM",
  "RF",
  "MAXNET",
  "SVM",
  "ENSEMBLE_MEAN",
  "ENSEMBLE_WEIGHTED_TSS"
)

for(alg in algorithms_to_plot) {
  
  # Africa figure
  fig_africa <- make_side_by_side_figure(
    algorithm_name = alg,
    region = "Africa",
    plot_fact = 1
  )
  
  if(!is.null(fig_africa)) {
    
    print(fig_africa)
    
    ggsave(
      file.path(fig_dir, paste0("FIGURE_", alg, "_TREATMENTS_AFRICA.png")),
      fig_africa,
      width = 13,
      height = 8,
      dpi = 400
    )
    
    ggsave(
      file.path(fig_dir, paste0("FIGURE_", alg, "_TREATMENTS_AFRICA.pdf")),
      fig_africa,
      width = 13,
      height = 8
    )
  }
  
  # Global figure
  fig_global <- make_side_by_side_figure(
    algorithm_name = alg,
    region = "Global",
    plot_fact = 2
  )
  
  if(!is.null(fig_global)) {
    
    print(fig_global)
    
    ggsave(
      file.path(fig_dir, paste0("FIGURE_", alg, "_TREATMENTS_GLOBAL.png")),
      fig_global,
      width = 13,
      height = 8,
      dpi = 400
    )
    
    ggsave(
      file.path(fig_dir, paste0("FIGURE_", alg, "_TREATMENTS_GLOBAL.pdf")),
      fig_global,
      width = 13,
      height = 8
    )
  }
}

message("Side-by-side treatment maps saved in: ", fig_dir)
# ============================================================
# 17. COUNTRY-LEVEL VALIDATION AGAINST KNOWN DISTRIBUTION
# Country-level Boyce-like index + zonal suitability
# ============================================================

library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(pROC)
library(rnaturalearth)

country_val_dir <- file.path(out_dir, "country_level_validation")
dir.create(country_val_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 1. Africa country polygons
# ------------------------------------------------------------

africa_sf <- rnaturalearth::ne_countries(
  continent = "Africa",
  returnclass = "sf",
  scale = "medium"
)

africa_sf <- st_make_valid(africa_sf)

# Standard country name field
africa_sf$Country <- africa_sf$admin

# ------------------------------------------------------------
# 2. Known occupied countries from published distribution map
# ------------------------------------------------------------

known_presence <- data.frame(
  Country_standard = c(
    "Angola",
    "Benin",
    "Botswana",
    "Cameroon",
    "Central African Republic",
    "Republic of the Congo",
    "Democratic Republic of the Congo",
    "Côte d'Ivoire",
    "Equatorial Guinea",
    "Ethiopia",
    "Gabon",
    "Ghana",
    "Guinea",
    "Guinea-Bissau",
    "Kenya",
    "Lesotho",
    "Liberia",
    "Mozambique",
    "Namibia",
    "Nigeria",
    "Senegal",
    "Sierra Leone",
    "Somalia",
    "South Africa",
    "Tanzania",
    "Togo",
    "Uganda",
    "Zambia",
    "Zimbabwe"
  ),
  Presence = 1
)
# ------------------------------------------------------------
# 2. Known occupied countries from published distribution map
# ------------------------------------------------------------

known_presence <- data.frame(
  Country_standard = c(
    "Angola", "Benin", "Botswana", "Cameroon",
    "Central African Republic", "Republic of the Congo",
    "Democratic Republic of the Congo", "Côte d'Ivoire",
    "Equatorial Guinea", "Ethiopia", "Gabon", "Ghana",
    "Guinea", "Guinea-Bissau", "Kenya", "Lesotho",
    "Liberia", "Mozambique", "Namibia", "Nigeria",
    "Senegal", "Sierra Leone", "Somalia", "South Africa",
    "Tanzania", "Togo", "Uganda", "Zambia", "Zimbabwe"
  ),
  Presence = 1
)

# Match names to Natural Earth admin names
africa_sf <- africa_sf %>%
  mutate(
    Country_standard = case_when(
      admin == "Ivory Coast" ~ "Côte d'Ivoire",
      admin == "Republic of the Congo" ~ "Republic of the Congo",
      admin == "Democratic Republic of the Congo" ~ "Democratic Republic of the Congo",
      admin == "United Republic of Tanzania" ~ "Tanzania",
      admin == "Central African Republic" ~ "Central African Republic",
      TRUE ~ admin
    )
  ) %>%
  left_join(known_presence, by = "Country_standard") %>%
  mutate(
    Presence = ifelse(is.na(Presence), 0, Presence)
  )

write.csv(
  st_drop_geometry(africa_sf) %>%
    dplyr::select(Country_standard, Presence),
  file.path(country_val_dir, "country_presence_absence_layer.csv"),
  row.names = FALSE
)
# ------------------------------------------------------------
# 3. Parse treatment and algorithm from filenames
# ------------------------------------------------------------

treatments_order <- c(
  "Random_background",
  "Kernel_background",
  "Persistence_zero_background",
  "Persistence_predictor",
  "Kernel_background_plus_predictor",
  "Persistence_zero_background_plus_predictor"
)

algorithms_all <- c(
  "GLM",
  "GAM",
  "RF",
  "MAXNET",
  "SVM",
  "ENSEMBLE_MEAN",
  "ENSEMBLE_WEIGHTED_TSS"
)

parse_map_info <- function(filename) {
  
  nm <- tools::file_path_sans_ext(basename(filename))
  nm <- gsub(paste0("^", species_name, "_"), "", nm)
  
  alg <- algorithms_all[sapply(algorithms_all, function(a) grepl(paste0("_", a, "$"), nm))]
  alg <- alg[which.max(nchar(alg))]
  
  tr <- gsub(paste0("_", alg, "$"), "", nm)
  
  data.frame(
    file = filename,
    Treatment = tr,
    Algorithm = alg
  )
}

map_files <- list.files(
  map_dir,
  pattern = "\\.tif$",
  full.names = TRUE
)

# Keep only maps produced during the treatment comparison
map_files <- map_files[
  !grepl("equal_contribution_ensemble", basename(map_files), ignore.case = TRUE)
]
map_info <- bind_rows(lapply(map_files, parse_map_info)) %>%
  filter(Treatment %in% treatments_order, Algorithm %in% algorithms_all)

# ------------------------------------------------------------
# 4. Zonal mean suitability per country
# ------------------------------------------------------------

africa_vect <- terra::vect(africa_sf)

country_stats_all <- list()

for(i in seq_len(nrow(map_info))) {
  
  f <- map_info$file[i]
  tr <- map_info$Treatment[i]
  alg <- map_info$Algorithm[i]
  
  message("Extracting country means: ", tr, " | ", alg)
  
  r <- terra::rast(f)
  
  # Crop to Africa for speed
  r_af <- terra::crop(r, africa_vect)
  
  z <- terra::extract(
    r_af,
    africa_vect,
    fun = mean,
    na.rm = TRUE
  )
  
  val_col <- names(z)[2]
  
  out <- st_drop_geometry(africa_sf) %>%
    mutate(
      Suitability_mean = z[[val_col]],
      Treatment = tr,
      Algorithm = alg
    ) %>%
    dplyr::select(
      Country = Country_standard,
      Presence,
      Treatment,
      Algorithm,
      Suitability_mean
    )
  
  country_stats_all[[paste(tr, alg, sep = "_")]] <- out
}

country_stats <- bind_rows(country_stats_all)

write.csv(
  country_stats,
  file.path(country_val_dir, "country_mean_suitability_all_models.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. Country-level Boyce-like index
# ------------------------------------------------------------

country_boyce <- function(suitability, presence, n_bins = 5) {
  
  df <- data.frame(
    suitability = suitability,
    presence = presence
  ) %>%
    filter(is.finite(suitability), !is.na(presence))
  
  if(length(unique(df$suitability)) < 3) return(NA_real_)
  if(sum(df$presence == 1) < 3) return(NA_real_)
  
  # Quantile bins avoid empty bins with small country-level sample
  probs <- seq(0, 1, length.out = n_bins + 1)
  breaks <- unique(as.numeric(quantile(df$suitability, probs = probs, na.rm = TRUE)))
  
  if(length(breaks) < 4) return(NA_real_)
  
  df$bin <- cut(
    df$suitability,
    breaks = breaks,
    include.lowest = TRUE,
    labels = FALSE
  )
  
  tab <- df %>%
    filter(!is.na(bin)) %>%
    group_by(bin) %>%
    summarise(
      midpoint = mean(suitability, na.rm = TRUE),
      expected = n() / nrow(df),
      observed = sum(presence == 1) / sum(df$presence == 1),
      .groups = "drop"
    ) %>%
    mutate(
      pe_ratio = observed / expected
    ) %>%
    filter(is.finite(pe_ratio), expected > 0)
  
  if(nrow(tab) < 3) return(NA_real_)
  
  suppressWarnings(
    cor(tab$midpoint, tab$pe_ratio, method = "spearman", use = "complete.obs")
  )
}

# ------------------------------------------------------------
# 6. Country-level validation summary
# ------------------------------------------------------------

country_validation_summary <- country_stats %>%
  group_by(Treatment, Algorithm) %>%
  summarise(
    Mean_suitability_occupied = mean(Suitability_mean[Presence == 1], na.rm = TRUE),
    Mean_suitability_background = mean(Suitability_mean[Presence == 0], na.rm = TRUE),
    Difference_occupied_minus_background =
      Mean_suitability_occupied - Mean_suitability_background,
    Ratio_occupied_to_background =
      Mean_suitability_occupied / Mean_suitability_background,
    Country_AUC = as.numeric(pROC::auc(Presence, Suitability_mean, quiet = TRUE)),
    Country_Boyce = country_boyce(Suitability_mean, Presence, n_bins = 5),
    Wilcoxon_p = suppressWarnings(
      wilcox.test(
        Suitability_mean[Presence == 1],
        Suitability_mean[Presence == 0]
      )$p.value
    ),
    N_occupied = sum(Presence == 1),
    N_background = sum(Presence == 0),
    .groups = "drop"
  ) %>%
  arrange(desc(Country_Boyce), desc(Difference_occupied_minus_background))

write.csv(
  country_validation_summary,
  file.path(country_val_dir, "country_validation_summary_boyce_like.csv"),
  row.names = FALSE
)

print(country_validation_summary)

# ------------------------------------------------------------
# 7. Treatment-level summary averaged across algorithms
#    INCLUDING WILCOXON SUMMARY
# ------------------------------------------------------------

treatment_country_summary <- country_validation_summary %>%
  group_by(Treatment) %>%
  summarise(
    Country_Boyce_mean = mean(Country_Boyce, na.rm = TRUE),
    Country_Boyce_sd   = sd(Country_Boyce, na.rm = TRUE),
    
    Country_AUC_mean = mean(Country_AUC, na.rm = TRUE),
    Country_AUC_sd   = sd(Country_AUC, na.rm = TRUE),
    
    Mean_occupied   = mean(Mean_suitability_occupied, na.rm = TRUE),
    Mean_background = mean(Mean_suitability_background, na.rm = TRUE),
    
    Difference_mean = mean(Difference_occupied_minus_background, na.rm = TRUE),
    Difference_sd   = sd(Difference_occupied_minus_background, na.rm = TRUE),
    
    Wilcoxon_p_median = median(Wilcoxon_p, na.rm = TRUE),
    Wilcoxon_p_min    = min(Wilcoxon_p, na.rm = TRUE),
    Wilcoxon_p_max    = max(Wilcoxon_p, na.rm = TRUE),
    
    Significant_models_p_lt_0.05 = sum(Wilcoxon_p < 0.05, na.rm = TRUE),
    N_models = n(),
    
    .groups = "drop"
  ) %>%
  arrange(desc(Country_Boyce_mean), desc(Difference_mean))

write.csv(
  treatment_country_summary,
  file.path(country_val_dir, "treatment_level_country_validation_summary_with_wilcoxon.csv"),
  row.names = FALSE
)

print(treatment_country_summary)

# ------------------------------------------------------------
# 8. Boxplot: occupied vs background countries
# ------------------------------------------------------------

p_country_box <- country_stats %>%
  mutate(
    Presence_class = ifelse(Presence == 1, "Known occupied countries", "Other African countries"),
    Treatment = factor(Treatment, levels = treatments_order)
  ) %>%
  ggplot(aes(x = Treatment, y = Suitability_mean, fill = Presence_class)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.85) +
  coord_flip() +
  theme_bw(base_size = 11) +
  labs(
    x = "Calibration strategy",
    y = "Country-level mean suitability",
    fill = NULL,
    title = "Country-level suitability in known occupied versus other African countries"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top"
  )

ggsave(
  file.path(country_val_dir, "Figure_country_suitability_boxplot.png"),
  p_country_box,
  width = 10,
  height = 6,
  dpi = 600
)

ggsave(
  file.path(country_val_dir, "Figure_country_suitability_boxplot.pdf"),
  p_country_box,
  width = 10,
  height = 6
)

# ------------------------------------------------------------
# 9. Barplot: Country-Boyce by treatment and algorithm
# ------------------------------------------------------------

p_boyce <- country_validation_summary %>%
  mutate(
    Treatment = factor(Treatment, levels = treatments_order)
  ) %>%
  ggplot(aes(x = Treatment, y = Country_Boyce, fill = Algorithm)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_bw(base_size = 11) +
  labs(
    x = "Calibration strategy",
    y = "Country-level Boyce-like index",
    title = "Country-level Boyce-like validation across treatments and algorithms"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

ggsave(
  file.path(country_val_dir, "Figure_country_boyce_by_treatment_algorithm.png"),
  p_boyce,
  width = 11,
  height = 7,
  dpi = 600
)

ggsave(
  file.path(country_val_dir, "Figure_country_boyce_by_treatment_algorithm.pdf"),
  p_boyce,
  width = 11,
  height = 7
)

message("Country-level validation completed. Outputs saved in: ", country_val_dir)

# ============================================================
# FIGURE: HEATMAP OF CROSS-VALIDATION PERFORMANCE
# TSS and AUC across calibration strategies and algorithms
# ============================================================

library(dplyr)
library(ggplot2)
library(patchwork)

# If not already loaded:
# summary_results <- read.csv(file.path(out_dir, "cross_validation_summary_with_ensembles.csv"))

heatmap_dir <- file.path(out_dir, "performance_heatmaps")
dir.create(heatmap_dir, showWarnings = FALSE, recursive = TRUE)

treatment_labels <- c(
  "Random_background" = "Random background",
  "Kernel_background" = "Kernel background",
  "Persistence_zero_background" = "Persistence background",
  "Persistence_predictor" = "Persistence predictor",
  "Kernel_background_plus_predictor" = "Kernel + predictor",
  "Persistence_zero_background_plus_predictor" = "Persistence background + predictor"
)

treatment_order <- c(
  "Persistence_zero_background_plus_predictor",
  "Persistence_zero_background",
  "Persistence_predictor",
  "Random_background",
  "Kernel_background_plus_predictor",
  "Kernel_background"
)

algorithm_order <- c(
  "GLM", "GAM", "RF", "MAXNET", "SVM",
  "ENSEMBLE_MEAN", "ENSEMBLE_WEIGHTED_TSS"
)

algorithm_labels <- c(
  "GLM" = "GLM",
  "GAM" = "GAM",
  "RF" = "RF",
  "MAXNET" = "MaxEnt",
  "SVM" = "SVM",
  "ENSEMBLE_MEAN" = "Mean ensemble",
  "ENSEMBLE_WEIGHTED_TSS" = "TSS-weighted ensemble"
)

plot_dat <- summary_results %>%
  mutate(
    Treatment_label = factor(
      treatment_labels[Treatment],
      levels = treatment_labels[treatment_order]
    ),
    Algorithm_label = factor(
      algorithm_labels[Algorithm],
      levels = algorithm_labels[algorithm_order]
    )
  )

# ----------------------------
# TSS heatmap
# ----------------------------

p_tss <- ggplot(
  plot_dat,
  aes(x = Algorithm_label, y = Treatment_label, fill = TSS_mean)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", TSS_mean)), size = 3) +
  scale_fill_viridis_c(
    name = "Mean TSS",
    option = "C",
    limits = c(0, 1)
  ) +
  theme_bw(base_size = 11) +
  labs(
    x = "Algorithm",
    y = "Calibration strategy",
    title = "Cross-validation TSS across calibration strategies and algorithms"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# ----------------------------
# AUC heatmap
# ----------------------------

p_auc <- ggplot(
  plot_dat,
  aes(x = Algorithm_label, y = Treatment_label, fill = AUC_mean)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", AUC_mean)), size = 3) +
  scale_fill_viridis_c(
    name = "Mean AUC",
    option = "C",
    limits = c(0.5, 1)
  ) +
  theme_bw(base_size = 11) +
  labs(
    x = "Algorithm",
    y = "Calibration strategy",
    title = "Cross-validation AUC across calibration strategies and algorithms"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# ----------------------------
# Combined figure
# ----------------------------

fig_heatmap <- p_tss / p_auc

ggsave(
  file.path(heatmap_dir, "Figure_cross_validation_TSS_AUC_heatmaps.png"),
  fig_heatmap,
  width = 12,
  height = 10,
  dpi = 600
)

ggsave(
  file.path(heatmap_dir, "Figure_cross_validation_TSS_AUC_heatmaps.pdf"),
  fig_heatmap,
  width = 12,
  height = 10
)

print(fig_heatmap)
# ============================================================
# FINAL PAPER FIGURE ONLY:
# MaxEnt + combined persistence-background and persistence predictor
# ============================================================

library(terra)
library(sf)
library(ggplot2)
library(tidyterra)
library(viridis)
library(rnaturalearth)
library(ggspatial)

# Install ggspatial if needed
# install.packages("ggspatial")

paper_fig_dir <- file.path(out_dir, "FINAL_PAPER_FIGURE")
dir.create(paper_fig_dir, showWarnings = FALSE, recursive = TRUE)

final_tif <- file.path(
  map_dir,
  paste0(
    species_name,
    "_Persistence_zero_background_plus_predictor_MAXNET.tif"
  )
)

if(!file.exists(final_tif)) {
  stop("Final MaxEnt raster not found: ", final_tif)
}

r_final <- terra::rast(final_tif)

world <- rnaturalearth::ne_countries(
  returnclass = "sf",
  scale = "medium"
)

# Optional reduce raster for faster plotting only
r_plot <- terra::aggregate(
  r_final,
  fact = 2,
  fun = mean,
  na.rm = TRUE
)

p_final_global <- ggplot() +
  geom_sf(
    data = world,
    fill = "grey96",
    color = NA
  ) +
  tidyterra::geom_spatraster(
    data = r_plot
  ) +
  geom_sf(
    data = world,
    fill = NA,
    color = "grey35",
    linewidth = 0.12
  ) +
  geom_point(
    data = occ,
    aes(x = lon, y = lat),
    size = 0.45,
    alpha = 0.75,
    color = "black"
  ) +
  scale_fill_viridis_c(
    name = "Suitability",
    option = "C",
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    labels = c("0", "0.25", "0.50", "0.75", "1.00"),
    na.value = "transparent",
    guide = guide_colorbar(
      barheight = unit(4.2, "cm"),
      barwidth = unit(0.35, "cm"),
      ticks = TRUE
    )
  ) +
  annotation_scale(
    location = "bl",
    width_hint = 0.22,
    text_cex = 0.65,
    line_width = 0.6
  ) +
  annotation_north_arrow(
    location = "bl",
    which_north = "true",
    pad_x = unit(0.15, "in"),
    pad_y = unit(0.55, "in"),
    style = north_arrow_fancy_orienteering,
    height = unit(0.75, "cm"),
    width = unit(0.75, "cm")
  ) +
  coord_sf(
    xlim = c(-180, 180),
    ylim = c(-60, 80),
    expand = FALSE
  ) +
  labs(
    title = expression(
      "Global climatic suitability of " * italic("Rhynchophorus phoenicis")
    ),
    subtitle = "Selected MaxEnt model calibrated with combined persistence-informed background and persistence predictor",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.text = element_text(color = "black", size = 8),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8),
    plot.margin = ggplot2::margin(8, 8, 8, 8)
  )

print(p_final_global)

ggsave(
  file.path(paper_fig_dir, "Figure_5_Final_MaxEnt_Global_Suitability.png"),
  p_final_global,
  width = 12,
  height = 6.5,
  dpi = 600
)

ggsave(
  file.path(paper_fig_dir, "Figure_5_Final_MaxEnt_Global_Suitability.pdf"),
  p_final_global,
  width = 12,
  height = 6.5
)

ggsave(
  file.path(paper_fig_dir, "Figure_5_Final_MaxEnt_Global_Suitability.tiff"),
  p_final_global,
  width = 12,
  height = 6.5,
  dpi = 600,
  compression = "lzw"
)

message("Final paper figure saved in: ", paper_fig_dir)
# ============================================================
# FINAL PAPER FIGURE ONLY:
# MaxEnt + combined persistence-background and persistence predictor
# No internal country borders; low suitability near white
# ============================================================

library(terra)
library(sf)
library(ggplot2)
library(tidyterra)
library(rnaturalearth)
library(ggspatial)

paper_fig_dir <- file.path(out_dir, "FINAL_PAPER_FIGURE")
dir.create(paper_fig_dir, showWarnings = FALSE, recursive = TRUE)

final_tif <- file.path(
  map_dir,
  paste0(species_name, "_Persistence_zero_background_plus_predictor_MAXNET.tif")
)

if(!file.exists(final_tif)) {
  stop("Final MaxEnt raster not found: ", final_tif)
}

r_final <- terra::rast(final_tif)

# Land outline only: removes internal country boundaries
land <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf"
)

land <- sf::st_union(sf::st_geometry(land))
land <- sf::st_sf(geometry = land)

# Reduce raster for plotting only
r_plot <- terra::aggregate(
  r_final,
  fact = 2,
  fun = mean,
  na.rm = TRUE
)

p_final_global <- ggplot() +
  
  tidyterra::geom_spatraster(data = r_plot) +
  
  geom_sf(
    data = land,
    fill = NA,
    colour = "grey35",
    linewidth = 0.22
  ) +
  
  geom_point(
    data = occ,
    aes(x = lon, y = lat),
    size = 0.35,
    shape = 21,
    fill = "black",
    colour = "white",
    stroke = 0.12,
    alpha = 0.8
  ) +
  
  scale_fill_stepsn(
    name = "Suitability",
    colours = c(
      "#ffffff",   # 0.00
      "#fff7ec",   # 0.05
      "#fee8c8",   # 0.10
      "#fdd49e",   # 0.20
      "#fdbb84",   # 0.40
      "#fc8d59",   # 0.60
      "#ef6548",   # 0.80
      "#d7301f",   # 0.90
      "#990000"    # 1.00
    ),
    breaks = c(
      0.00,
      0.05,
      0.10,
      0.20,
      0.40,
      0.60,
      0.80,
      0.90,
      1.00
    ),
    limits = c(0,1),
    na.value = "white",
    show.limits = TRUE,
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(8, "cm"),
      barheight = unit(0.35, "cm"),
      ticks = TRUE
    )
  ) +
  
  annotation_scale(
    location = "bl",
    width_hint = 0.22,
    text_cex = 0.65,
    line_width = 0.6,
    pad_x = unit(0.20, "in"),
    pad_y = unit(0.15, "in")
  ) +
  
  annotation_north_arrow(
    location = "bl",
    which_north = "true",
    pad_x = unit(0.18, "in"),
    pad_y = unit(0.55, "in"),
    style = north_arrow_fancy_orienteering,
    height = unit(0.75, "cm"),
    width = unit(0.75, "cm")
  ) +
  
  coord_sf(
    xlim = c(-180, 180),
    ylim = c(-60, 80),
    expand = FALSE
  ) +
  
  
  theme_bw(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.15),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.text = element_text(color = "black", size = 8),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 5),
    plot.margin = ggplot2::margin(8, 8, 8, 8)
  )

print(p_final_global)

ggsave(
  file.path(paper_fig_dir, "Figure_5_Final_MaxEnt_Global_Suitability.png"),
  p_final_global,
  width = 12,
  height = 6.5,
  dpi = 600
)

ggsave(
  file.path(paper_fig_dir, "Figure_5_Final_MaxEnt_Global_Suitability.pdf"),
  p_final_global,
  width = 12,
  height = 6.5
)

ggsave(
  file.path(paper_fig_dir, "Figure_5_Final_MaxEnt_Global_Suitability.tiff"),
  p_final_global,
  width = 12,
  height = 6.5,
  dpi = 600,
  compression = "lzw"
)

message("Final paper figure saved in: ", paper_fig_dir)

