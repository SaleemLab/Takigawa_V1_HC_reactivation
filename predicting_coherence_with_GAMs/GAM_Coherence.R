# ==============================================================================
# V1-HC COHERENCE: STAGE 5 - FINAL MODEL & DIAGNOSTICS
# ==============================================================================

# --- 1. Load Required Libraries ---
packages <- c("mgcv", "gratia", "ggplot2", "dplyr", "tidyr")
for (p in packages) {
  if (!require(p, character.only = TRUE)) install.packages(p)
}

library(mgcv)
library(gratia)
library(ggplot2)
library(dplyr)
library(tidyr)

# --- 2. Load and Prepare Data ---
message("Loading data...")

### Load -200 to 200ms or Load -100 to 100ms
dat <- read.csv("C:/Users/masah/Desktop/V1_HC_reactivation/predicting_coherence_with_GAMs/v1_hc_data_geo_-200_200ms.csv")
dat <- read.csv("C:/Users/masah/Desktop/V1_HC_reactivation/predicting_coherence_with_GAMs/v1_hc_data_geo_-100_100ms.csv")

dat$SessionID <- as.factor(dat$SessionID)

dat$HPC_logodds_Z <- as.numeric(scale(dat$HPC_logodds_Z))
dat$V1_logodds_Z <- as.numeric(scale(dat$V1_logodds_Z))
dat$V1_logodds_PRE_Z <- as.numeric(scale(dat$V1_logodds_PRE_Z))

# --- 3. Clean Data ---
z_thresh <- 2.56  # include <99th centiles of data

message(sprintf("Trimming extreme outliers beyond +/- %s Z-scores...", z_thresh))
dat_clean <- dat %>%
  filter(
    RipplePower_Z > -z_thresh & RipplePower_Z < z_thresh,
    SpindlePower_Match_Z > -z_thresh & SpindlePower_Match_Z < z_thresh,
    SpindlePower_NonMatch_Z > -z_thresh & SpindlePower_NonMatch_Z < z_thresh,
  )
  # Standardize the existing Coherence column from the CSV

# Standardize the Dependent Variable (Coherence)
# This makes 'Partial Effect' units equal to Standard Deviations of Coherence
  #mutate(Event_Coherence_Post_Z = as.numeric(scale(Event_Coherence_Post)))

# ==============================================================================
# --- TASK 5.1: FIT THE MINIMAL ADEQUATE MODEL ---
# ==============================================================================
message("\nFitting the Final Minimal Adequate Model...")


mdl_final <- bam(Event_Coherence_Post_GeoMean ~ 
                   s(RipplePower_Z, k = 5) + 
                   s(SpindlePower_Match_Z, k = 5) + 
                   s(SpindlePower_NonMatch_Z, k = 5) +
                   #te(SOPhase_Match, SOPhase_NonMatch, bs = "cc", k = 8) + 
                   ti(SOPhase_Match, SOPhase_NonMatch, bs = "cc", k = 8) + 
                   s(SOPhase_Match, bs = "cc", k = 8) + 
                   s(SOPhase_NonMatch, bs = "cc", k = 8) + 

                   # 4. Control Term
                   s(AnimalID, bs = "re") +
                   s(SessionID, bs = "re"), 
                 
                 data = dat_clean, 
                 method = "fREML", 
                 discrete = TRUE, 
                 nthreads = 4)

message("\n--- FINAL MODEL SUMMARY ---")
print(summary(mdl_final))


# ==============================================================================
# --- TASK 5.2: FORMAL DIAGNOSTICS (GAM.CHECK) ---
# ==============================================================================
message("\n--- RUNNING DIAGNOSTICS (gam.check) ---")
# gam.check(mdl_final)

# ==============================================================================
# --- STAGE 5 VISUALIZATIONS: PUBLICATION FIGURES ---
# ==============================================================================
message("\nGenerating Final Publication Figures...")

# 1. Define the folder path you want to go to
my_folder <- "C:/Users/masah/Desktop/V1_HC_reactivation/predicting_coherence_with_GAMs/Model1"
#my_folder <- "C:/Users/masah/Desktop/V1_HC_reactivation/predicting_coherence_with_GAMs/final_model"
#my_folder <- "C:/Users/masah/Desktop/V1_HC_reactivation/predicting_coherence_with_GAMs/final_model_200ms"


# 2. Check if it exists; if not, create it
if (!dir.exists(my_folder)) {
  dir.create(my_folder, recursive = TRUE)
}

# 3. Change the working directory to that folder
setwd(my_folder)


# 1. Main Effect: Ripple Power
# cairo_pdf("ripple_power_effect.pdf", width = 4.3, height = 4.3)
# p_ripple <- draw(mdl_final, select = "s(RipplePower_Z)", residuals = FALSE, rug = FALSE) + 
#   theme_bw(base_family = "Arial") + 
#   theme(aspect.ratio = 1) +
#   labs(title = "Main Effect: Ripple Power", x = "Ripple Power (Z)", y = "Partial Effect")
# 
# print(p_ripple)
# dev.off()

# Calculate scaling factors
raw_breaks <- c(5, 10, 15)
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean$RipplePower_Z[which.min(abs(dat_clean$RipplePower - val))]
})

cairo_pdf("ripple_power_effect_raw.pdf", width = 4.3, height = 4.3)
p_ripple_raw <- draw(mdl_final, select = "s(RipplePower_Z)", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  labs(
    title = "Ripple Power and coherence", 
    x = "Ripple Power (Raw Units)", 
    y = "Partial Effect"
  )

print(p_ripple_raw)
dev.off()



# 2. Main Effect: Spindle Power

# cairo_pdf("spindle_power_matched_effect.pdf", width = 4.3, height = 4.3)
# p_spindle <- draw(mdl_final, select = "s(SpindlePower_Match_Z)", residuals = FALSE, rug = FALSE) + 
#   theme_bw(base_family = "Arial") + 
#   theme(aspect.ratio = 1) +
#   labs(title = "Main Effect: Spindle Power Matched", x = "Spindle Power Matched (Z)", y = "Partial Effect")
# 
# print(p_spindle)
# dev.off()

######## Matched (Dominant hemisphere) Spindle Power
# Calculate scaling factors
raw_breaks_spindle <- c(-1, 1, 3)
z_breaks_spindle <- sapply(raw_breaks_spindle, function(val) {
  dat_clean$SpindlePower_Match_Z[which.min(abs(dat_clean$SpindlePower_Match - val))]
})

cairo_pdf("spindle_power_matched_effect_raw.pdf", width = 4.3, height = 4.3)
p_spindle_raw <- draw(mdl_final, select = "s(SpindlePower_Match_Z)", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks_spindle, labels = raw_breaks_spindle) +
  labs(
    title = "Spindle Power Matched and coherence", 
    x = "Spindle Power (Raw Units)", 
    y = "Partial Effect"
  )

print(p_spindle_raw)
dev.off()


######## Non-matched (non-Dominant hemisphere) spindle power

cairo_pdf("spindle_power_Nonmatched_effect.pdf", width = 4.3, height = 4.3)
p_spindle <- draw(mdl_final, select = "s(SpindlePower_Match_Z)", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  labs(title = "Main Effect: Spindle Power Matched", x = "Spindle Power Matched (Z)", y = "Partial Effect")

print(p_spindle)
dev.off()

# Calculate scaling factors
raw_breaks_spindle <- c(-1, 1, 3)
z_breaks_spindle <- sapply(raw_breaks_spindle, function(val) {
  dat_clean$SpindlePower_NonMatch_Z[which.min(abs(dat_clean$SpindlePower_NonMatch - val))]
})

cairo_pdf("spindle_power_Nonmatched_effect_raw.pdf", width = 4.3, height = 4.3)
p_spindle_raw <- draw(mdl_final, select = "s(SpindlePower_NonMatch_Z)", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks_spindle, labels = raw_breaks_spindle) +
  labs(
    title = "Spindle Power NonMatched and coherence", 
    x = "Spindle Power (Raw Units)", 
    y = "Partial Effect"
  )

print(p_spindle_raw)
dev.off()



# ---------------------------------------------------------
# Plot 1: Main Effect of Match Phase
# ---------------------------------------------------------
sm_match <- smooth_estimates(mdl_final, smooth = "s(SOPhase_Match)", n = 200) %>%
  mutate(
    # The wrap-around trick!
    SOPhase_Match = ifelse(SOPhase_Match < 0, SOPhase_Match + 2*pi, SOPhase_Match),
    .lower_ci = .estimate - (1.96 * .se),
    .upper_ci = .estimate + (1.96 * .se)
  ) %>% arrange(SOPhase_Match)

# dev.new(noRStudioGD = TRUE)
p_match <- ggplot(sm_match, aes(x = SOPhase_Match, y = .estimate)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, fill = "dodgerblue") +
  geom_line(linewidth = 1, color = "dodgerblue") +
  theme_bw() + 
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi)),
                     limits = c(0, 2*pi), expand = c(0, 0)) +
  labs(title = "Isolated Main Effect: Match SO Phase", 
       x = "Match SO Phase (Rad)", y = "Partial Effect")

cairo_pdf("SO_Match_effect.pdf", width = 4.3, height = 4.3)
print(p_match)
dev.off()
# ---------------------------------------------------------
# Plot 2: Main Effect of Non-Match Phase
# ---------------------------------------------------------
sm_nonmatch <- smooth_estimates(mdl_final, smooth = "s(SOPhase_NonMatch)", n = 200) %>%
  mutate(
    # The wrap-around trick!
    SOPhase_NonMatch = ifelse(SOPhase_NonMatch < 0, SOPhase_NonMatch + 2*pi, SOPhase_NonMatch),
    .lower_ci = .estimate - (1.96 * .se),
    .upper_ci = .estimate + (1.96 * .se)
  ) %>% arrange(SOPhase_NonMatch)

# dev.new(noRStudioGD = TRUE)
p_nonmatch <- ggplot(sm_nonmatch, aes(x = SOPhase_NonMatch, y = .estimate)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, fill = "firebrick") +
  geom_line(linewidth = 1, color = "firebrick") +
  theme_bw() + 
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi)),
                     limits = c(0, 2*pi), expand = c(0, 0)) +
  labs(title = "Isolated Main Effect: Non-Match SO Phase", 
       x = "Non-Match SO Phase (Rad)", y = "Partial Effect")
print(p_nonmatch)

cairo_pdf("SO_NonMatch_effect.pdf", width = 4.3, height = 4.3)
print(p_nonmatch)
dev.off()

# ---------------------------------------------------------
# Plot 3: The "Pure" Phase Interaction (ti term)
# ---------------------------------------------------------

# Note: Ensure your model (mdl_final) uses these renamed terms
sm_2d <- smooth_estimates(mdl_final, smooth = "ti(SOPhase_Match,SOPhase_NonMatch)", n = 100) %>%
  mutate(
    .lower_ci = .estimate - (1.96 * .se),
    .upper_ci = .estimate + (1.96 * .se)
  )

p_synergy <- sm_2d %>%
  mutate(
    SOPhase_Match = ifelse(SOPhase_Match < 0, SOPhase_Match + 2*pi, SOPhase_Match),
    SOPhase_NonMatch = ifelse(SOPhase_NonMatch < 0, SOPhase_NonMatch + 2*pi, SOPhase_NonMatch)
  ) %>%
  arrange(SOPhase_Match, SOPhase_NonMatch) %>%
  ggplot(aes(x = SOPhase_Match, y = SOPhase_NonMatch, fill = .estimate)) +
  geom_tile(width = (2*pi)/99, height = (2*pi)/99) + 
  geom_contour(aes(z = .estimate), color = "black", alpha = 0.2) + 
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", midpoint = 0, name = "Effect") +
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  scale_y_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  theme_minimal(base_family = "Arial") +
  theme(aspect.ratio = 1) +
  labs(title = "Inter-hemispheric Phase interaction effect on coherence",
       x = "Match SO Phase (Rad)", y = "Non-Match SO Phase (Rad)")

# dev.new(noRStudioGD = TRUE)
cairo_pdf("SO_interaction_effect.pdf", width = 5.2, height = 4.3)
print(p_synergy)
dev.off()

# ==============================================================================
# --- RECREATING THE TOTAL LANDSCAPE ---
# ==============================================================================

# 1. Create a clean 2D grid
pred_grid_phase <- expand.grid(
  SOPhase_Match = seq(0, 2*pi, length.out = 100),
  SOPhase_NonMatch = seq(0, 2*pi, length.out = 100),
  RipplePower_Z = 0,
  SpindlePower_Match_Z = 0,
  SpindlePower_NonMatch_Z = 0,
  # Fixed the ID assignment
  SessionID = dat_clean$SessionID[1],
  AnimalID = dat_clean$AnimalID[1] 
)

# 2. Extract terms 
term_preds_phase <- predict(mdl_final, newdata = pred_grid_phase, type = "terms")

# 3. Sum the specific terms
# Note: Ensure these strings match your mdl_final formula exactly
pred_grid_phase$Reconstructed_TE <- 
  term_preds_phase[, "s(SOPhase_Match)"] + 
  term_preds_phase[, "s(SOPhase_NonMatch)"] + 
  term_preds_phase[, "ti(SOPhase_Match,SOPhase_NonMatch)"]

message("\nReconstructing Total Phase Synergy from Decomposed Terms...")

p_reconstructed <- pred_grid_phase %>%
  mutate(
    SOPhase_Match = ifelse(SOPhase_Match < 0, SOPhase_Match + 2*pi, SOPhase_Match),
    SOPhase_NonMatch = ifelse(SOPhase_NonMatch < 0, SOPhase_NonMatch + 2*pi, SOPhase_NonMatch)
  ) %>%
  arrange(SOPhase_Match, SOPhase_NonMatch) %>%
  ggplot(aes(x = SOPhase_Match, y = SOPhase_NonMatch, fill = Reconstructed_TE)) +
  geom_tile(width = (2*pi)/99, height = (2*pi)/99) + 
  geom_contour(aes(z = Reconstructed_TE), color = "black", alpha = 0.2) + 
  scale_fill_gradient2(
    low = "dodgerblue", 
    mid = "white", 
    high = "firebrick", 
    midpoint = 0, 
    name = "Effect"
  ) +
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  scale_y_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  theme_minimal(base_family = "Arial") +
  theme(aspect.ratio = 1) +
  labs(
    title = "Inter-hemispheric Phase effect on coherence", 
    subtitle = "Reconstructed: s(Match) + s(NonMatch) + ti(Match, NonMatch)",
    x = "Match SO Phase (Rad)", 
    y = "Non-Match SO Phase (Rad)"
  )

# dev.new(noRStudioGD = TRUE)
cairo_pdf("SO_total_effect.pdf", width = 5.2, height = 4.3)
print(p_reconstructed)
dev.off()


# ==============================================================================
# --- 3. MULTI-METRIC EFFECT SIZE CALCULATIONS ---
# ==============================================================================
message("\nCalculating 4 Effect Size Metrics (This may take a minute)...")
# 
# # EXACT formula components to rebuild the models robustly
# formula_terms <- c(
#   "s(RipplePower_Z, k = 5)", 
#   "s(SpindlePower_Match_Z, k = 5)", 
#   "s(SpindlePower_NonMatch_Z, k = 5)", 
#   "ti(SpindlePower_Match_Z,SpindlePower_NonMatch_Z, k = 5)",
#   "s(SOPhase_Match, bs = 'cc', k = 8)",
#   "s(SOPhase_NonMatch, bs = 'cc', k = 8)",
#   "ti(SOPhase_Match,SOPhase_NonMatch, bs = 'cc', k = 8)",
#   "ti(RipplePower_Z,SOPhase_Match, bs = c('tp', 'cc'), k = c(5, 8))"
# )
# 
# 
# # EXACT labels output by summary() and smooth_estimates()
# smooth_labels <- c(
#   "s(RipplePower_Z)", 
#   "s(SpindlePower_Match_Z)", 
#   "s(SpindlePower_NonMatch_Z)",
#   "ti(SpindlePower_Match_Z,SpindlePower_NonMatch_Z)",
#   "s(SOPhase_Match)",
#   "s(SOPhase_NonMatch)",
#   "ti(SOPhase_Match,SOPhase_NonMatch)",
#   "ti(RipplePower_Z,SOPhase_Match)"
# )


# EXACT formula components to rebuild the models robustly
formula_terms <- c(
  "s(RipplePower_Z, k = 5)",
  "s(SpindlePower_Match_Z, k = 5)",
  "s(SpindlePower_NonMatch_Z, k = 5)",
  "s(SOPhase_Match, bs = 'cc', k = 8)",
  "s(SOPhase_NonMatch, bs = 'cc', k = 8)",
  "ti(SOPhase_Match,SOPhase_NonMatch, bs = 'cc', k = 8)"
)


# EXACT labels output by summary() and smooth_estimates()
smooth_labels <- c(
  "s(RipplePower_Z)",
  "s(SpindlePower_Match_Z)",
  "s(SpindlePower_NonMatch_Z)",
  "s(SOPhase_Match)",
  "s(SOPhase_NonMatch)",
  "ti(SOPhase_Match,SOPhase_NonMatch)"
)

library(parallel)
library(pbapply)
library(dplyr)
library(stringr)
library(tidyr)

B <- 1000  # Set to 1000 for your final analysis run
RE_TERMS <- c("s(SessionID, bs = 're')", "s(AnimalID, bs = 're')")

message(sprintf("\nLaunching %d Case Bootstrap Replicates...", B))

run_one_bootstrap <- function(rep_id, original_data, formula_terms, smooth_labels) {
  
  # 1. Resample data with replacement
  boot_data <- original_data[sample(nrow(original_data), replace = TRUE), ]
  
  # 2. Fit Full Model
  full_form <- as.formula(paste(
    "Event_Coherence_Post_GeoMean ~",
    paste(c(formula_terms, RE_TERMS), collapse = " + ")
  ))
  
  mdl_f <- mgcv::bam(full_form, data = boot_data, method = "fREML", discrete = TRUE)
  
  f_sum <- summary(mdl_f)
  full_dev <- f_sum$dev.expl * 100
  res_df <- mdl_f$df.residual
  sum_tab <- as.data.frame(f_sum$s.table)
  sum_tab$Term <- rownames(sum_tab)
  
  # Pre-allocate container for this run's metrics
  run_res <- data.frame(Term = smooth_labels, Partial_Deviance = NA, Eta_Sq_Partial = NA, Amplitude = NA, RMS = NA, Rep = rep_id)
  
  # 3. Process every smooth term
  for(j in 1:length(smooth_labels)) {
    # A. Partial Deviance (Drop-One refit)
    act_terms <- c(formula_terms[-j], "s(SessionID, bs = 're')")
    r_form <- as.formula(paste("Event_Coherence_Post_GeoMean ~", paste(act_terms, collapse = " + ")))
    mdl_r <- mgcv::bam(r_form, data = boot_data, method = "fREML", discrete = TRUE)
    run_res$Partial_Deviance[j] <- full_dev - (summary(mdl_r)$dev.expl * 100)
    
    # B. Partial Eta-Squared
    t_row <- sum_tab[sum_tab$Term == smooth_labels[j], ]
    if(nrow(t_row) == 1) {
      run_res$Eta_Sq_Partial[j] <- (t_row$F * t_row$edf) / ((t_row$F * t_row$edf) + res_df)
    }
    
    # C & D. Amplitude and RMS via completely self-contained custom grid
    vars_in_smooth <- all.vars(as.formula(paste("~", smooth_labels[j])))
    vars_in_smooth <- vars_in_smooth[!vars_in_smooth %in% c("s", "ti", "bs", "k")]
    
    # Build a clean evaluation sequence for variables in this smooth
    slice_args <- lapply(vars_in_smooth, function(v) {
      seq(min(boot_data[[v]], na.rm = TRUE), max(boot_data[[v]], na.rm = TRUE), length.out = 50)
    })
    names(slice_args) <- vars_in_smooth
    
    # Use base R expand.grid to keep data structure strictly local
    grid_clean <- do.call(expand.grid, slice_args)
    
    # Add dummy/mean columns for any other variables that predict() might structurally look for
    all_model_vars <- all.vars(full_form)[-1] 
    missing_vars <- setdiff(all_model_vars, colnames(grid_clean))
    
    for(mv in missing_vars) {
      if(mv == "SessionID") {
        grid_clean[[mv]] <- boot_data$SessionID[1] 
      } else if(is.numeric(boot_data[[mv]])) {
        grid_clean[[mv]] <- mean(boot_data[[mv]], na.rm = TRUE)
      } else {
        grid_clean[[mv]] <- boot_data[[mv]][1]
      }
    }
    
    # Isolate target smooth prediction using the lpmatrix
    Xp <- predict(mdl_f, newdata = grid_clean, type = "lpmatrix")
    smooth_cols <- grep(smooth_labels[j], colnames(Xp), fixed = TRUE)
    
    Xp_isolated <- matrix(0, nrow = nrow(Xp), ncol = ncol(Xp))
    Xp_isolated[, smooth_cols] <- Xp[, smooth_cols]
    
    fit_isolated <- Xp_isolated %*% coef(mdl_f)
    run_res$Amplitude[j] <- max(fit_isolated) - min(fit_isolated)
    run_res$RMS[j]       <- sqrt(mean(fit_isolated^2))
  }
  
  return(run_res)
}

# Run loop with progress bar
set.seed(42)
boot_results_list <- pblapply(1:B, run_one_bootstrap, 
                              original_data = dat_clean, 
                              formula_terms = formula_terms, 
                              smooth_labels = smooth_labels)

# Clean out any NULLs if any structural issues happen
boot_results_list <- boot_results_list[!sapply(boot_results_list, is.null)]

# Combine results safely
boot_df <- do.call(rbind, boot_results_list)

# 4. Compute Median and 95% CIs from the Bootstrap Distribution
final_dashboard_data <- boot_df %>%
  dplyr::group_by(Term) %>%
  dplyr::summarise(
    dplyr::across(c(Partial_Deviance, Eta_Sq_Partial, Amplitude, RMS),
                  list(Val = ~median(.x, na.rm = TRUE),
                       Lwr = ~quantile(.x, probs = 0.025, na.rm = TRUE),
                       Upr = ~quantile(.x, probs = 0.975, na.rm = TRUE)),
                  .names = "{.col}__{.fn}")
  )

# ==============================================================================
# --- 3E. SAVE RAW BOOTSTRAP DISTRIBUTIONS (ALL REPLICATES) ---
# ==============================================================================
message("Saving raw bootstrap iterations for distribution archives...")

# Format the column headers cleanly before exporting
raw_iterations_clean <- boot_df %>%
  select(Rep, Term, Partial_Deviance, Eta_Sq_Partial, Amplitude, RMS) %>%
  rename(
    Replicate        = Rep,
    Deviance_Value   = Partial_Deviance,
    Eta_Sq_Value     = Eta_Sq_Partial,
    Amplitude_Value  = Amplitude,
    RMS_Value        = RMS
  ) %>%
  arrange(Term, Replicate)

# Export the raw simulation matrix to a CSV file
write.csv(raw_iterations_clean, "GAM_model_raw_bootstrap_iterations.csv", row.names = FALSE)

message("Raw bootstrap iteration file saved successfully: 'GAM_model_raw_bootstrap_iterations.csv'")


library(tidyverse)

# ==============================================================================
# --- 4. GRAPH BOOTSTRAPPED EFFECT SIZE DASHBOARD (BAR PLOT EDITION) ---
# ==============================================================================

# 1. Pivot the Median Estimates
eff_val <- final_dashboard_data %>%
  select(Term, ends_with("__Val")) %>%
  rename_with(~str_remove(., "__Val"), -Term) %>%
  pivot_longer(cols = -Term, names_to = "Metric", values_to = "Value")

# 2. Pivot the Lower Confidence Bounds
eff_lwr <- final_dashboard_data %>%
  select(Term, ends_with("__Lwr")) %>%
  rename_with(~str_remove(., "__Lwr"), -Term) %>%
  pivot_longer(cols = -Term, names_to = "Metric", values_to = "Lower")

# 3. Pivot the Upper Confidence Bounds
eff_upr <- final_dashboard_data %>%
  select(Term, ends_with("__Upr")) %>%
  rename_with(~str_remove(., "__Upr"), -Term) %>%
  pivot_longer(cols = -Term, names_to = "Metric", values_to = "Upper")

# 4. Merge everything cleanly and apply factor labels
plot_data <- eff_val %>%
  left_join(eff_lwr, by = c("Term", "Metric")) %>%
  left_join(eff_upr, by = c("Term", "Metric")) %>%
  mutate(Metric = factor(Metric, 
                         levels = c("Partial_Deviance", "Eta_Sq_Partial", "Amplitude", "RMS"),
                         labels = c("Deviance Explained (%)", "Partial Eta-Squared", 
                                    "Peak-to-Trough Amplitude", "RMS Effect")))

# 5. Construct the Custom Bar Plot
p_bars_with_ci <- ggplot(plot_data, aes(x = reorder(Term, Value), y = Value, fill = Term)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", alpha = 0.5) +
  geom_bar(stat = "identity", show.legend = FALSE, alpha = 0.85, width = 0.75) +
  # Asymmetric error bars marking the lower and upper bounds cleanly
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.25, color = "black", size = 0.6) +
  coord_flip() +
  facet_wrap(~Metric, scales = "free_x", ncol = 2) +
  scale_fill_viridis_d(option = "mako", direction = -1) +
  theme_bw() +
  labs(title = "GAMM Effect Size & Variance Metrics",
       subtitle = "Bars denote bootstrap medians; error bars represent 95% non-parametric CIs.",
       x = NULL, y = NULL) +
  theme(
    strip.text = element_text(face = "bold", size = 11, colour = "black"),
    strip.background = element_rect(fill = "gray95", color = "gray80"),
    axis.text.y = element_text(face = "plain", size = 10, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray95"),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    panel.spacing = unit(1.2, "lines")
  )

# Render to RStudio console
print(p_bars_with_ci)

# 6. Save vector graphic output
cairo_pdf("Model_Effect_Sizes_With_CI.pdf", width = 9.5, height = 6.5)
print(p_bars_with_ci)
dev.off()



# ==============================================================================
# --- 5. SAVE COMPREHENSIVE MODEL OUTPUT WITH 95% CIs ---
# ==============================================================================
library(dplyr)
library(tidyr)
library(stringr)

# 1. Format the basic model summary parameters
model_stats <- as.data.frame(summary(mdl_final)$s.table) %>%
  mutate(Term = rownames(.))

# 2. Flatten the nested stats into individual columns
flat_bootstrap_results <- final_dashboard_data %>%
  pivot_longer(cols = -Term, names_to = "Combined", values_to = "Value") %>%
  separate(Combined, into = c("Metric", "Stat"), sep = "__") %>%
  mutate(New_Col_Name = paste0(Metric, "_", Stat)) %>%
  select(-Metric, -Stat) %>%
  pivot_wider(names_from = New_Col_Name, values_from = Value)

# 3. Join the parametric summary table with the resampled metrics
combined_results <- model_stats %>%
  left_join(flat_bootstrap_results, by = "Term") %>%
  select(
    Term, edf, Ref.df, F, `p-value`,
    Partial_Deviance_Val, Partial_Deviance_Lwr, Partial_Deviance_Upr,
    Eta_Sq_Partial_Val, Eta_Sq_Partial_Lwr, Eta_Sq_Partial_Upr,
    Amplitude_Val, Amplitude_Lwr, Amplitude_Upr,
    RMS_Val, RMS_Lwr, RMS_Upr
  ) %>%
  rename(
    p_value              = `p-value`,
    Deviance_Median      = Partial_Deviance_Val,
    Deviance_CI_Lower    = Partial_Deviance_Lwr,
    Deviance_CI_Upper    = Partial_Deviance_Upr,
    Eta_Sq_Median        = Eta_Sq_Partial_Val,
    Eta_Sq_CI_Lower      = Eta_Sq_Partial_Lwr,
    Eta_Sq_CI_Upper      = Eta_Sq_Partial_Upr,
    Amplitude_Median     = Amplitude_Val,
    Amplitude_CI_Lower   = Amplitude_Lwr,
    Amplitude_CI_Upper   = Amplitude_Upr,
    RMS_Median           = RMS_Val,
    RMS_CI_Lower         = RMS_Lwr,
    RMS_CI_Upper         = RMS_Upr
  )

# 4. Export the comprehensive data table
write.csv(combined_results, "GAM_model_CI_output.csv", row.names = FALSE)

message("Master spreadsheet file saved successfully: 'GAM_model_CI_output.csv'")
message("\nAnalysis Pipeline Complete!")
