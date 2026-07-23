# V1-HC Coordinated Reactivation Code

This repository contains the analysis scripts to replicate the core findings from:

> **Interplay of sleep neural oscillations enhances coordinated memory reactivation between cortex and hippocampus**
> *Masahiro Takigawa, D. Tong, E. A. B. Horrocks, A. B. Saleem, and D. Bendor*
> bioRxiv, 2026. DOI: [10.64898/2026.06.12.731367v1](https://www.biorxiv.org/content/10.64898/2026.06.12.731367v1)

## Code strcture

### Visualisation of track-selective activity in V1 and HC during VR track running
*   [V1_HC_Track_FR_distribution.m](file:///C:/Users/masah/Desktop/V1_HC_reactivation/V1_HC_Track_FR_distribution.m): Calculation and visualisation of spatial map and firing rate distribution on Track L and Track R for left and right V1 and hippocampus (HPC)
*   [PLS_ROC_RUN_analysis.m](file:///C:/Users/masah/Desktop/V1_HC_reactivation/PLS_ROC_RUN_analysis.m): Plots ROC curves (AUC) representing track discrimination performance during RUN sessions using Partial Least Squares (PLS) decoders.

### Visualisation of spatial PSTH and peri-ripple PSTH
*   [interactive_V1_ripple_explorer.m]: A MATLAB App designed to interactively explore V1 neural firing profiles, spatial tuning curves, and ripple PSTHs relative to track differences and LME results.

### Reactivation Coherence & Sleep Oscillation AUC Traces
These scripts calculate a sliding-window Area Under the Curve (AUC) metric to show how much V1's reactivation bias aligns with HPC's bias under various oscillatory conditions:
*   [reactivation_AUC_coherence_all_events.m]: Baseline temporal AUC for all sleep ripple events.
*   [reactivation_AUC_coherence_ripple_power.m]: Evaluates V1-HPC coherence across different ripple power percentiles.
*   [reactivation_AUC_coherence_spindle_power.m]: Evaluates coherence across different cortical spindle power bands (dominant/matched vs. non-matched hemispheres).
*   [reactivation_AUC_coherence_SO_phase.m](file:///C:/Users/masah/Desktop/V1_HC_reactivation/reactivation_AUC_coherence_SO_phase.m): Evaluates coherence during slow oscillation peaks vs. troughs.
*   [reactivation_AUC_coherence_SO_trough_uni_vs_bilateral.m]: Compares V1-HPC reactivation coherence under unilateral vs. bilateral slow oscillation troughs.

### GAMM Coherence Modeling (`/predicting_coherence_with_GAMs`)
*   [GAM_Coherence.R]: Fits Generalized Additive Mixed Models (GAMMs) using `mgcv::bam` to model non-linear oscillatory effects (SO phase, spindle power, ripple power) on V1-HC reactivation cohernece (Geometric Mean of V1 and HC track bias) with random effects for Animal ID and Session ID. Includes a non-parametric case bootstrapping module to yield 95% confidence intervals for effect sizes.


## Instruction
All processed data needed to run the code is found within `/processed_data` folder. Simply download this github repository to run the code for visualization.

## Expected Run Time

- **Plotting scripts:** typically under 30 seconds.
- **Scripts involving shuffling** (e.g., permutation tests): roughly 10–30 minutes, depending on data size and number of iterations.

## Requirements
* Tested using MATLAB 2024a (https://uk.mathworks.com/products/matlab.html)
* Tested using RStudio 2026.04.0 (https://posit.co/download/rstudio-desktop/)