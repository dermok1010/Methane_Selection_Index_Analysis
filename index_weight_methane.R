# ============================================================
# SELECTION RESPONSE COMPARISON OF METHANE TRAITS
# (point estimates only - Monte Carlo removed)
# Uses ASReml .pvc estimates
# ============================================================
setwd("/home/dermot.kelly/Dermot_analysis/Phd/Paper_3/")

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(ggtext)
})

# -----------------------------
# 0) Names of traits + means for ratio Taylor
# -----------------------------
X <- "ch4_g_day2_1v3"
Y <- "liveweight"
Z <- "adg"

mu_X <- 17.90444
mu_Y <- 62.72

# -----------------------------
# 1) Point estimates of (co)variance components from .pvc outputs
# -----------------------------

# ---- Bivariate: ch4_g_day2_1v3 (trait1) vs liveweight (trait2) ----
# From bi_ch4_weight.pvc
VA_X_hat   <- 2.44586     # Additive genetic variance: CH4
VA_Y_hat   <- 57.7836     # Additive genetic variance: liveweight
rg_XY_hat  <- 0.5018      # Genetic correlation: CH4-liveweight

VE_X_hat   <- 10.6096     # Residual variance: CH4
VE_Y_hat   <- 14.8424     # Residual variance: liveweight
re_XY_hat  <- 0.2397      # Residual correlation: CH4-liveweight

VPE_XY_hat <- 2.66683     # PE variance: ide(ANI_ID), CH4-liveweight bivariate

# ---- Bivariate: liveweight (trait1) vs adg_g (trait2) ----
# From bi_weight_adg.pvc
rg_YZ_hat  <- 0.7763      # Genetic correlation: liveweight-ADG

# ---- Bivariate: ch4_g_day2_1v3 (trait1) vs adg_g (trait2) ----
# From bi_ch4_adg_g.pvc
VA_Z_from_X_hat <- 2420.85   # Additive genetic variance: ADG (g/day)
VE_Z_from_X_hat <- 134.237   # Residual variance: ADG (g/day)
rg_XZ_hat       <- 0.4048    # Genetic correlation: CH4-ADG
VPE_XZ_hat      <- 1.96238   # PE variance: ide(ANI_ID), CH4-ADG bivariate

##########
# Lambda / coefficients of additive genetic variation
##########
CVA_X  <- sqrt(VA_X_hat) / mu_X
CVA_Y  <- sqrt(VA_Y_hat) / mu_Y
lambda <- CVA_X / CVA_Y

CVA_X
CVA_Y
lambda

# -----------------------------
# 2) Helpers
# -----------------------------
# Returns the smallest eigenvalue of a matrix.
eig_min <- function(M) min(eigen(M, symmetric = TRUE, only.values = TRUE)$values)

# Checks if a matrix is valid (positive definite) and repairs it if not.
make_PD <- function(M){
  if(eig_min(M) < 1e-10) as.matrix(nearPD(M)$mat) else M
}

# Reconstructs covariance from variances and correlation.
cov_from_r <- function(Vi, Vj, r) r * sqrt(Vi * Vj)

# Computes selection index response for all grid directions simultaneously.
calc_resp_2trait <- function(G2, P2, A){
  B     <- solve(P2, G2 %*% A)
  denom <- sqrt(colSums(B * (P2 %*% B)))
  R     <- (t(G2) %*% B) / matrix(denom, nrow = 2, ncol = length(denom), byrow = TRUE)
  list(B = B, denom = denom, R = R)
}

# Computes correlated genetic response in traits not in the index (e.g. ADG).
calc_resp_H_from_B <- function(G_HI, B, denom){
  (G_HI %*% B) / matrix(denom, nrow = nrow(G_HI), ncol = length(denom), byrow = TRUE)
}


get_hull <- function(df, x = "lw", y = "ch4"){
  if(nrow(df) < 3) return(df[0, , drop = FALSE])
  idx <- chull(df[[x]], df[[y]])
  df[idx, , drop = FALSE]
}

# -----------------------------
# 3) Grid of economic weights
# -----------------------------
w_CH4_seq <- seq(-10, 10, by = 0.5)
w_LW_seq  <- seq(-10, 10, by = 0.5)

grid <- expand.grid(w_CH4 = w_CH4_seq, w_LW = w_LW_seq) %>%
  filter(!(w_CH4 == 0 & w_LW == 0))

A_grid <- rbind(grid$w_CH4, grid$w_LW)
rownames(A_grid) <- c(X, Y)

# ── SPECIAL POINT 1: RATIO (TAYLOR LINEARISATION) ────────────────────────────
# Minimise CH4/liveweight. Taylor expansion gradient at trait means:
#   a_CH4 = -1/mu_Y
#   a_LW  = +mu_X/mu_Y^2
a_ratio <- matrix(c(-1/mu_Y, mu_X/(mu_Y^2)), ncol = 1)
rownames(a_ratio) <- c(X, Y)

# -----------------------------
# 4) Baseline matrices (CH4-liveweight bivariate)
# -----------------------------
CovA_XY_hat <- cov_from_r(VA_X_hat, VA_Y_hat, rg_XY_hat)
CovE_XY_hat <- cov_from_r(VE_X_hat, VE_Y_hat, re_XY_hat)
CovP_XY_hat <- CovA_XY_hat + CovE_XY_hat    # no PE covariance

VP_X_hat <- VA_X_hat + VE_X_hat + VPE_XY_hat
VP_Y_hat <- VA_Y_hat + VE_Y_hat + VPE_XY_hat

G2_hat <- matrix(c(VA_X_hat,    CovA_XY_hat,
                   CovA_XY_hat, VA_Y_hat),
                 nrow = 2, byrow = TRUE, dimnames = list(c(X,Y), c(X,Y)))

P2_hat <- matrix(c(VP_X_hat,    CovP_XY_hat,
                   CovP_XY_hat, VP_Y_hat),
                 nrow = 2, byrow = TRUE, dimnames = list(c(X,Y), c(X,Y)))

G2_hat <- make_PD(G2_hat)
P2_hat <- make_PD(P2_hat)

out_hat <- calc_resp_2trait(G2_hat, P2_hat, A_grid)

resp_all <- grid %>%
  mutate(ch4 = as.numeric(out_hat$R[1,]),
         lw  = as.numeric(out_hat$R[2,]))

hull_df        <- get_hull(resp_all, x = "lw", y = "ch4") %>% mutate(type = "Frontier")
hull_df_closed <- bind_rows(hull_df, hull_df[1, , drop = FALSE])

# Ratio special point
out_ratio_hat <- calc_resp_2trait(G2_hat, P2_hat, a_ratio)
ratio_point <- tibble(type = "Ratio (Taylor)",
                      ch4  = as.numeric(out_ratio_hat$R[1,1]),
                      lw   = as.numeric(out_ratio_hat$R[2,1]))

# Residual special point (CH4 adjusted for liveweight)
beta_P_hat  <- as.numeric(P2_hat[X,Y] / P2_hat[Y,Y])
a_res_hat   <- matrix(c(-1, beta_P_hat), ncol = 1); rownames(a_res_hat) <- c(X,Y)
out_res_hat <- calc_resp_2trait(G2_hat, P2_hat, a_res_hat)
resid_point <- tibble(type = paste0("Residual (beta_P=", signif(beta_P_hat, 4), ")"),
                      ch4  = as.numeric(out_res_hat$R[1,1]),
                      lw   = as.numeric(out_res_hat$R[2,1]))

# ── 5) BUILD THE 3x3 G MATRIX FOR ADG ────────────────────────────────────────
# Diagonal sources:
#   CH4 (X): variance from the CH4-liveweight bivariate
#   liveweight (Y): variance from the CH4-liveweight bivariate
#   ADG (Z): variance from the CH4-ADG bivariate

VA_Z3_hat <- VA_Z_from_X_hat

CovA_XZ_hat <- cov_from_r(VA_X_hat, VA_Z3_hat, rg_XZ_hat)
CovA_YZ_hat <- cov_from_r(VA_Y_hat, VA_Z3_hat, rg_YZ_hat)

G3_hat <- matrix(0, 3, 3, dimnames = list(c(X,Y,Z), c(X,Y,Z)))
G3_hat[X,X] <- VA_X_hat
G3_hat[Y,Y] <- VA_Y_hat
G3_hat[Z,Z] <- VA_Z3_hat
G3_hat[X,Y] <- CovA_XY_hat; G3_hat[Y,X] <- CovA_XY_hat
G3_hat[X,Z] <- CovA_XZ_hat; G3_hat[Z,X] <- CovA_XZ_hat
G3_hat[Y,Z] <- CovA_YZ_hat; G3_hat[Z,Y] <- CovA_YZ_hat
G3_hat <- make_PD(G3_hat)

G_HI_hat <- G3_hat[, c(X,Y), drop = FALSE]
RH_hat   <- calc_resp_H_from_B(G_HI_hat, out_hat$B, out_hat$denom)

resp_all_adg <- resp_all %>% mutate(adg = as.numeric(RH_hat[3,]))

lr_hat <- resp_all_adg %>% filter(lw > 0, ch4 < 0)

# ADG at the two special points
adg_ratio_hat <- as.numeric(calc_resp_H_from_B(G_HI_hat, out_ratio_hat$B, out_ratio_hat$denom)[3,1])
adg_resid_hat <- as.numeric(calc_resp_H_from_B(G_HI_hat, out_res_hat$B,   out_res_hat$denom)[3,1])

special_adg <- tibble(
  point        = c("Ratio (Taylor)", paste0("Residual (beta_P=", signif(beta_P_hat, 4), ")")),
  ch4_response = c(ratio_point$ch4, resid_point$ch4),
  lw_response  = c(ratio_point$lw,  resid_point$lw),
  adg_response = c(adg_ratio_hat,   adg_resid_hat)
)

# -----------------------------
# 6) Diagnostic plots
# -----------------------------
ggplot() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(data = resp_all, aes(x = lw, y = ch4), alpha = 0.12, size = 1) +
  geom_path(data = hull_df_closed, aes(x = lw, y = ch4), linewidth = 1.0) +
  geom_point(data = ratio_point, aes(x = lw, y = ch4), size = 4, shape = 17) +
  geom_point(data = resid_point, aes(x = lw, y = ch4), size = 4, shape = 15) +
  labs(
    x        = "Predicted response in liveweight (i=1 normalised)",
    y        = "Predicted response in CH4/day (i=1 normalised)",
    title    = "2D response region and frontier for CH4 & liveweight",
    subtitle = "Line: convex-hull frontier. Triangle: ratio (Taylor). Square: residual (X-betaY)."
  ) +
  theme_minimal()

# -----------------------------
# 7) Results summary (point estimates)
# -----------------------------
special_points <- special_adg %>%
  transmute(point, ch4 = ch4_response, lw = lw_response, adg = adg_response)

cat("\n=============================\n")
cat("CORE COMPARISON: Ratio vs Residual (point estimates)\n")
cat("=============================\n")
print(special_points)

adg_range_lr <- lr_hat %>%
  summarise(
    n       = n(),
    adg_min = min(adg, na.rm = TRUE),
    adg_max = max(adg, na.rm = TRUE),
    adg_p05 = quantile(adg, 0.05, na.rm = TRUE),
    adg_p50 = quantile(adg, 0.50, na.rm = TRUE),
    adg_p95 = quantile(adg, 0.95, na.rm = TRUE)
  )

cat("\n=============================\n")
cat("ADG RANGE UNDER CH4-LW INDEX (favourable quadrant)\n")
cat("=============================\n")
print(adg_range_lr)

lr_growth_ok <- lr_hat %>% filter(adg >= 0)

cat("\n=============================\n")
cat("GROWTH-PROTECTED REGION (ADG >= 0 within favourable quadrant)\n")
cat("=============================\n")
cat("Count:", nrow(lr_growth_ok), "out of", nrow(lr_hat), "\n")

if(nrow(lr_growth_ok) > 0){
  cat("\nMax methane reduction achievable with ADG >= 0 (within grid):\n")
  print(lr_growth_ok %>% arrange(ch4) %>% slice(1) %>% select(ch4, lw, adg, w_CH4, w_LW))
  
  cat("\nMax ADG achievable while still in favourable quadrant and ADG >= 0:\n")
  print(lr_growth_ok %>% arrange(desc(adg)) %>% slice(1) %>% select(ch4, lw, adg, w_CH4, w_LW))
} else {
  cat("No growth-protected points found on this grid.\n")
}


# ============================================================
# 8) FRONTIER PLOT — HULL LINE COLOURED BY ADG RESPONSE
# ============================================================

lr <- resp_all_adg

hull_lr_adg <- resp_all_adg %>%
  mutate(angle = atan2(ch4, lw),
         r     = sqrt(lw^2 + ch4^2)) %>%
  group_by(cut(angle, 720)) %>%          # 720 bins = 0.5 degree resolution
  slice_max(r, n = 1) %>%
  ungroup() %>%
  arrange(angle)

# Fill any floating point misses via nearest neighbour
if (any(is.na(hull_lr_adg$adg))) {
  missing_idx <- which(is.na(hull_lr_adg$adg))
  for (i in missing_idx) {
    dists <- sqrt(
      (resp_all_adg$lw  - hull_lr_adg$lw[i])^2 +
        (resp_all_adg$ch4 - hull_lr_adg$ch4[i])^2
    )
    hull_lr_adg$adg[i] <- resp_all_adg$adg[which.min(dists)]
  }
}

adg_lims     <- range(lr$adg, na.rm = TRUE)
pct_positive <- round(100 * mean(lr$adg >= 0, na.rm = TRUE), 1)
med_adg      <- median(lr$adg, na.rm = TRUE)

palette_stops <- c("#8B0000", "#CC3300", "#FF6644", "#FAEBD7",
                   "#88CC88", "#2E8B2E", "#145214")

value_stops <- rescale(c(
  adg_lims[1],
  adg_lims[1] * 0.5,
  adg_lims[1] * 0.1,
  0,
  adg_lims[2] * 0.1,
  adg_lims[2] * 0.5,
  adg_lims[2]
))

adg_col_scale <- scale_colour_gradientn(
  colours = palette_stops,
  values  = value_stops,
  limits  = adg_lims,
  breaks  = c(adg_lims[1], adg_lims[2]),
  labels  = sprintf("%.2f", c(adg_lims[1], adg_lims[2])),
  name    = expression("ADG correlated response (g day"^{-1}*")")
)

p_frontier <- ggplot() +
  
  # Layer 1: Background cloud
  geom_point(
    data  = lr,
    aes(x = lw, y = ch4, colour = adg),
    alpha = 0.25,
    size  = 1.4
  ) +
  adg_col_scale +
  
  # Layer 2: Hull line coloured by ADG
  geom_path(
    data      = hull_lr_adg,
    aes(x = lw, y = ch4, colour = adg),
    linewidth = 2.2,
    lineend   = "round",
    linejoin  = "round"
  ) +
  
  # Layer 3: Ratio special point
  geom_point(
    data   = ratio_point,
    aes(x = lw, y = ch4),
    size   = 5,
    shape  = 24,
    fill   = "white",
    colour = "#1A1A2E",
    stroke = 1.5
  ) +
  annotate(
    "richtext",
    x     = ratio_point$lw - 0.25,
    y     = ratio_point$ch4 - 0.3,
    label = paste0(
      "<b>Ratio</b><br>",
      "\u0394CH4: ", round(ratio_point$ch4, 2), "<br>",
      "\u0394LW: ",  round(ratio_point$lw,  2), "<br>",
      "\u0394ADG: ", round(adg_ratio_hat,    2)
    ),
    size        = 4.3,
    colour      = "#1A1A2E",
    hjust       = 0,
    fill        = NA,
    label.color = NA
  ) +
  
  # Layer 4: Residual special point
  geom_point(
    data   = resid_point %>% filter(lw > 0, ch4 < 0),
    aes(x = lw, y = ch4),
    size   = 5,
    shape  = 22,
    fill   = "white",
    colour = "#1A1A2E",
    stroke = 1.5
  ) +
  annotate(
    "richtext",
    x     = resid_point$lw + 0.3,
    y     = resid_point$ch4 - 0.2,
    label = paste0(
      "<b>Residual</b><br>",
      "\u0394CH4: ", round(resid_point$ch4, 2), "<br>",
      "\u0394LW: ",  round(resid_point$lw,  2), "<br>",
      "\u0394ADG: ", round(adg_resid_hat,    2)
    ),
    size        = 4.3,
    colour      = "#1A1A2E",
    hjust       = 0,
    fill        = NA,
    label.color = NA
  ) +
  
  # Reference lines
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.35) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.35) +
  
  labs(
    x = expression("Live weight (kg)"),
    y = expression("CH"[4]*" production (g day"^{-1}*")")
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_blank(),
    plot.subtitle    = element_blank(),
    axis.title       = element_text(size = 10, colour = "#1A1A2E"),
    legend.position  = "bottom",
    legend.title     = element_text(size = 9),
    legend.key.width = unit(2.5, "cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey93")
  ) +
  coord_cartesian(xlim = c(0, NA))

p_frontier

ggsave(
  "selection_index_P3/outputs/frontier_adg_liveweight.png",
  plot   = p_frontier,
  width  = 12,
  height = 8,
  dpi    = 600,
  bg     = "white"
)

# ------------------------------------------------------------
# Land-use variant (wider x range, no special points annotated)
# ------------------------------------------------------------
land_use_frontier <- ggplot() +
  geom_point(
    data  = lr,
    aes(x = lw, y = ch4, colour = adg),
    alpha = 0.25,
    size  = 1.4
  ) +
  adg_col_scale +
  geom_path(
    data      = hull_lr_adg,
    aes(x = lw, y = ch4, colour = adg),
    linewidth = 2.2,
    lineend   = "round",
    linejoin  = "round"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.35) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.35) +
  labs(
    x = expression("Live weight (kg)"),
    y = expression("CH"[4]*" production (g day"^{-1}*")")
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_blank(),
    plot.subtitle    = element_blank(),
    axis.title       = element_text(size = 10, colour = "#1A1A2E"),
    legend.position  = "bottom",
    legend.title     = element_text(size = 9),
    legend.key.width = unit(2.5, "cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey93")
  ) +
  coord_cartesian(xlim = c(-2, NA))

land_use_frontier

ggsave(
  "selection_index_P3/outputs/land_use.png",
  plot   = land_use_frontier,
  width  = 12,
  height = 8,
  dpi    = 600,
  bg     = "white"
)