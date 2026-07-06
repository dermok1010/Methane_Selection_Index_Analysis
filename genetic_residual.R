# ============================================================
# SELECTION RESPONSE COMPARISON OF METHANE TRAITS WITH MONTE CARLO UNCERTAINTY
# Uses ASReml .pvc estimates (estimate followed by SE)
#
# MODIFICATION: Added genetic residual index (beta_G = CovA(X,Y)/VA(Y))
# alongside the existing phenotypic residual (beta_P = CovP(X,Y)/VP(Y)).
# The genetic residual targets animals with lower CH4 than expected given
# their GENETIC merit for MBW, removing the genetic covariance between
# CH4 and MBW rather than the phenotypic association.
# ============================================================
setwd("/home/dermot.kelly/Dermot_analysis/Phd/Paper_3/")
asreml_data <- read.csv("genetic_analysis/asreml_scripts/P3_data.csv")

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

packageVersion("Matrix")
packageVersion("dplyr")
packageVersion("tidyr")
packageVersion("ggplot2")

citation("Matrix")
citation("dplyr")
citation("tidyr")
citation("ggplot2")

# -----------------------------
# 0) Names of traits + means for ratio Taylor
# -----------------------------
X <- "ch4_g_day2_1v3"
Y <- "Metabolic_BW"
Z <- "adg"

mu_X <- 17.90444
mu_Y <- 21.96824

# -----------------------------
# 1) Point estimates of (co)variance components + SEs from .pvc outputs
# -----------------------------
# ---- Bivariate: Metabolic_BW (trait1) vs ch4_g_day2_1v3 (trait2) ----
# Additive
VA_Y_hat  <- 1.51857   ; SE_VA_Y  <- 0.164883    # Additive genetic variance: MBW
VA_X_hat  <- 2.72527   ; SE_VA_X  <- 0.251641    # Additive genetic variance: CH4
rg_XY_hat <- -0.033    ; SE_rg_XY <- 0.0838      # Genetic correlation: MBW-CH4

# Residual
VE_Y_hat  <- 0.964220  ; SE_VE_Y  <- 0.168305E-01  # Residual variance: MBW
VE_X_hat  <- 10.7364   ; SE_VE_X  <- 0.170771      # Residual variance: CH4
re_XY_hat <- 0.2358    ; SE_re_XY <- 0.0115         # Residual correlation: MBW-CH4

# PE: ide(ANI_ID) shared across both traits
VPE_XY_hat <- 2.38080  ; SE_VPE_XY <- 0.157983     # Permanent environment variance

# ---- Bivariate: Metabolic_BW (trait1) vs adg (trait2) ----
VA_Y2_hat  <- 4.29765        ; SE_VA_Y2  <- 0.890151E-01
VA_Z_from_Y_hat <- 0.194251E-02 ; SE_VA_Z_from_Y <- 0.194836E-03
rg_YZ_hat  <- 0.6878         ; SE_rg_YZ  <- 0.0341
VE_Y2_hat  <- 0.994608       ; SE_VE_Y2  <- 0.174554E-01
VE_Z_from_Y_hat <- 0.134158E-03 ; SE_VE_Z_from_Y <- 0.516191E-05
re_YZ_hat  <- 0.0830         ; SE_re_YZ  <- 0.0343
VPE_YZ_hat <- 0.793797E-03   ; SE_VPE_YZ <- 0.167115E-03

# ---- Bivariate: ch4_g_day2_1v3 (trait1) vs adg (trait2) ----
VA_X2_hat  <- 4.62718        ; SE_VA_X2  <- 0.230782
VA_Z_from_X_hat <- 0.149254E-02 ; SE_VA_Z_from_X <- 0.205584E-03
rg_XZ_hat  <- 0.4104         ; SE_rg_XZ  <- 0.0486
VE_X2_hat  <- 10.9767        ; SE_VE_X2  <- 0.176163
VE_Z_from_X_hat <- 0.133262E-03 ; SE_VE_Z_from_X <- 0.510974E-05
re_XZ_hat  <- 0.0846         ; SE_re_XZ  <- 0.0300
VPE_XZ_hat <- 0.817222E-03   ; SE_VPE_XZ <- 0.183645E-03

# -----------------------------
# 2) Helpers
# -----------------------------
eig_min <- function(M) min(eigen(M, symmetric=TRUE, only.values=TRUE)$values)

make_PD <- function(M){
  if(eig_min(M) < 1e-10) as.matrix(nearPD(M)$mat) else M
}

rpos_lognormal <- function(n, x, se){
  if(is.na(se) || se <= 0) return(rep(x, n))
  cv <- se / x
  sigma <- sqrt(log(1 + cv^2))
  mu <- log(x) - 0.5 * sigma^2
  exp(rnorm(n, mean=mu, sd=sigma))
}

rcor_fisher <- function(n, r, se_r){
  if(is.na(se_r) || se_r <= 0) return(rep(r, n))
  r <- max(min(r, 0.999), -0.999)
  z <- atanh(r)
  se_z <- se_r / (1 - r^2)
  z_draw <- rnorm(n, mean=z, sd=se_z)
  pmax(pmin(tanh(z_draw), 0.999), -0.999)
}

cov_from_r <- function(Vi, Vj, r){
  r * sqrt(Vi * Vj)
}

compute_b <- function(G, P, a) solve(P, G %*% a)

calc_resp_2trait <- function(G2, P2, A){
  B <- solve(P2, G2 %*% A)
  denom <- sqrt(colSums(B * (P2 %*% B)))
  R <- (t(G2) %*% B) / matrix(denom, nrow=2, ncol=length(denom), byrow=TRUE)
  list(B=B, denom=denom, R=R)
}

calc_resp_H_from_B <- function(G_HI, B, denom){
  (G_HI %*% B) / matrix(denom, nrow=nrow(G_HI), ncol=length(denom), byrow=TRUE)
}

get_hull <- function(df, x="mbw", y="ch4"){
  if(nrow(df) < 3) return(df[0, , drop=FALSE])
  idx <- chull(df[[x]], df[[y]])
  df[idx, , drop=FALSE]
}

mc_summ <- function(x){
  x <- x[is.finite(x)]
  tibble(
    mean = mean(x),
    mc_sd = sd(x),
    lwr_2.5 = as.numeric(quantile(x, 0.025)),
    upr_97.5 = as.numeric(quantile(x, 0.975)),
    n = length(x)
  )
}

# -----------------------------
# 3) Grid (fixed)
# -----------------------------
w_CH4_seq <- seq(-10, 10, by=0.5)
w_MBW_seq <- seq(-10, 10, by=0.5)

grid <- expand.grid(w_CH4 = w_CH4_seq, w_MBW = w_MBW_seq) %>%
  filter(!(w_CH4 == 0 & w_MBW == 0))

A_grid <- rbind(grid$w_CH4, grid$w_MBW)
rownames(A_grid) <- c(X, Y)

# ── SPECIAL POINT 1: RATIO (TAYLOR LINEARISATION) ────────────────────────────
a_ratio <- matrix(c(-1/mu_Y, mu_X/(mu_Y^2)), ncol=1)
rownames(a_ratio) <- c(X, Y)

# -----------------------------
# 4) Baseline matrices
# -----------------------------
CovA_XY_hat <- cov_from_r(VA_X_hat, VA_Y_hat, rg_XY_hat)
CovE_XY_hat <- cov_from_r(VE_X_hat, VE_Y_hat, re_XY_hat)
CovP_XY_hat <- CovA_XY_hat + CovE_XY_hat
VP_X_hat <- VA_X_hat + VE_X_hat + VPE_XY_hat
VP_Y_hat <- VA_Y_hat + VE_Y_hat + VPE_XY_hat

G2_hat <- matrix(c(VA_X_hat, CovA_XY_hat,
                   CovA_XY_hat, VA_Y_hat),
                 nrow=2, byrow=TRUE, dimnames=list(c(X,Y), c(X,Y)))

P2_hat <- matrix(c(VP_X_hat, CovP_XY_hat,
                   CovP_XY_hat, VP_Y_hat),
                 nrow=2, byrow=TRUE, dimnames=list(c(X,Y), c(X,Y)))

G2_hat <- make_PD(G2_hat)
P2_hat <- make_PD(P2_hat)

out_hat <- calc_resp_2trait(G2_hat, P2_hat, A_grid)

resp_all <- grid %>%
  mutate(ch4 = as.numeric(out_hat$R[1,]),
         mbw = as.numeric(out_hat$R[2,]))

hull_df <- get_hull(resp_all, x="mbw", y="ch4") %>% mutate(type="Frontier")
hull_df_closed <- bind_rows(hull_df, hull_df[1, , drop=FALSE])

# ── RATIO POINT ───────────────────────────────────────────────────────────────
out_ratio_hat <- calc_resp_2trait(G2_hat, P2_hat, a_ratio)
ratio_point <- tibble(type="Ratio (Taylor)",
                      ch4 = as.numeric(out_ratio_hat$R[1,1]),
                      mbw = as.numeric(out_ratio_hat$R[2,1]))

# ── PHENOTYPIC RESIDUAL POINT ─────────────────────────────────────────────────
# Selects on CH4 after removing its phenotypic regression on MBW.
# beta_P = CovP(X,Y) / VP(Y)
# Captures animals with lower CH4 than expected for their observed body size.
beta_P_hat <- as.numeric(P2_hat[X,Y] / P2_hat[Y,Y])
a_res_hat <- matrix(c(-1, beta_P_hat), ncol=1); rownames(a_res_hat) <- c(X,Y)
out_res_hat <- calc_resp_2trait(G2_hat, P2_hat, a_res_hat)
resid_point <- tibble(type=paste0("Residual-P (beta_P=", signif(beta_P_hat,4), ")"),
                      ch4 = as.numeric(out_res_hat$R[1,1]),
                      mbw = as.numeric(out_res_hat$R[2,1]))

# ── GENETIC RESIDUAL POINT ────────────────────────────────────────────────────
# NEW: Selects on CH4 after removing its GENETIC regression on MBW.
# beta_G = CovA(X,Y) / VA(Y)
#
# Interpretation: this index targets animals whose breeding value for CH4
# is lower than expected given their breeding value for MBW. It removes the
# genetic covariance between CH4 and MBW rather than the total phenotypic
# association. This is conceptually "purer" in that it operates on the
# latent genetic architecture, but it still uses phenotypic measurements
# in the index weights (via P^{-1}).
#
# The desired direction vector is a = [-1, beta_G]:
#   - penalise CH4 breeding value
#   - give credit for MBW breeding value in proportion to how much MBW
#     genetically predicts CH4
# Note: beta_G and beta_P will differ whenever the genetic and residual
# correlations between CH4 and MBW differ — which they do here (rg = -0.033
# vs re = 0.236), so the two residual strategies target meaningfully
# different aspects of the genetic architecture.
beta_G_hat <- as.numeric(G2_hat[X,Y] / G2_hat[Y,Y])
a_gres_hat <- matrix(c(-1, beta_G_hat), ncol=1); rownames(a_gres_hat) <- c(X,Y)
out_gres_hat <- calc_resp_2trait(G2_hat, P2_hat, a_gres_hat)
gresid_point <- tibble(type=paste0("Residual-G (beta_G=", signif(beta_G_hat,4), ")"),
                       ch4 = as.numeric(out_gres_hat$R[1,1]),
                       mbw = as.numeric(out_gres_hat$R[2,1]))

# ── SECTION 5: BUILD THE 3x3 G MATRIX FOR ADG ────────────────────────────────
VA_X3_hat <- VA_X_hat
VE_X3_hat <- VE_X_hat
VPE_X3_hat <- VPE_XY_hat

VA_Y3_hat <- VA_Y_hat
VE_Y3_hat <- VE_Y_hat
VPE_Y3_hat <- VPE_XY_hat

VA_Z3_hat <- VA_Z_from_X_hat
VE_Z3_hat <- VE_Z_from_X_hat
VPE_Z3_hat <- VPE_XZ_hat

CovA_XZ_hat <- cov_from_r(VA_X3_hat, VA_Z3_hat, rg_XZ_hat)
CovA_YZ_hat <- cov_from_r(VA_Y3_hat, VA_Z3_hat, rg_YZ_hat)

G3_hat <- matrix(0, 3, 3, dimnames=list(c(X,Y,Z), c(X,Y,Z)))
G3_hat[X,X] <- VA_X3_hat
G3_hat[Y,Y] <- VA_Y3_hat
G3_hat[Z,Z] <- VA_Z3_hat
G3_hat[X,Y] <- CovA_XY_hat; G3_hat[Y,X] <- CovA_XY_hat
G3_hat[X,Z] <- CovA_XZ_hat; G3_hat[Z,X] <- CovA_XZ_hat
G3_hat[Y,Z] <- CovA_YZ_hat; G3_hat[Z,Y] <- CovA_YZ_hat
G3_hat <- make_PD(G3_hat)

G_HI_hat <- G3_hat[, c(X,Y), drop=FALSE]
RH_hat <- calc_resp_H_from_B(G_HI_hat, out_hat$B, out_hat$denom)

resp_all_adg <- resp_all %>% mutate(adg = as.numeric(RH_hat[3,]) * 1000)

lr_hat <- resp_all_adg %>% filter(mbw > 0, ch4 < 0)

# -----------------------------
# 6) Baseline plots
# -----------------------------
ggplot() +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  geom_point(data=resp_all, aes(x=mbw, y=ch4), alpha=0.12, size=1) +
  geom_path(data=hull_df_closed, aes(x=mbw, y=ch4), linewidth=1.0) +
  geom_point(data=ratio_point, aes(x=mbw, y=ch4), size=4, shape=17) +
  geom_point(data=resid_point, aes(x=mbw, y=ch4), size=4, shape=15) +
  geom_point(data=gresid_point, aes(x=mbw, y=ch4), size=4, shape=23) +
  labs(
    x="Predicted response in Metabolic_BW (i=1 normalised)",
    y="Predicted response in CH4/day (i=1 normalised)",
    title="2D response region and frontier for CH4 & Metabolic_BW",
    subtitle="Line: frontier. Triangle: ratio. Square: phenotypic residual. Diamond: genetic residual."
  ) +
  theme_minimal()

ggplot() +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  geom_point(data=resp_all, aes(x=mbw, y=ch4), alpha=0.12, size=1) +
  geom_path(data=hull_df_closed, aes(x=mbw, y=ch4), linewidth=1.0) +
  geom_point(data=ratio_point, aes(x=mbw, y=ch4), size=4, shape=17) +
  geom_point(data=resid_point, aes(x=mbw, y=ch4), size=4, shape=15) +
  geom_point(data=gresid_point, aes(x=mbw, y=ch4), size=4, shape=23) +
  coord_cartesian(
    xlim = c(0, max(resp_all$mbw, na.rm=TRUE)),
    ylim = c(min(resp_all$ch4, na.rm=TRUE), 0)
  ) +
  labs(
    x="Predicted response in Metabolic_BW (i=1 normalised)",
    y="Predicted response in CH4/day (i=1 normalised)",
    title="2D response region and frontier for CH4 & Metabolic_BW",
    subtitle="Zoomed to bottom-right quadrant (↑MBW, ↓CH4). Diamond: genetic residual."
  ) +
  theme_minimal()

# -----------------------------
# 7) Monte Carlo
# -----------------------------
set.seed(1)
n_sim <- 50

draw_VA_X  <- rpos_lognormal(n_sim, VA_X_hat,  SE_VA_X)
draw_VA_Y  <- rpos_lognormal(n_sim, VA_Y_hat,  SE_VA_Y)
draw_VE_X  <- rpos_lognormal(n_sim, VE_X_hat,  SE_VE_X)
draw_VE_Y  <- rpos_lognormal(n_sim, VE_Y_hat,  SE_VE_Y)
draw_VPE_XY <- rpos_lognormal(n_sim, VPE_XY_hat, SE_VPE_XY)

draw_rg_XY <- rcor_fisher(n_sim, rg_XY_hat, SE_rg_XY)
draw_re_XY <- rcor_fisher(n_sim, re_XY_hat, SE_re_XY)

draw_VA_Z <- rpos_lognormal(n_sim, VA_Z_from_X_hat, SE_VA_Z_from_X)
draw_VE_Z <- rpos_lognormal(n_sim, VE_Z_from_X_hat, SE_VE_Z_from_X)
draw_VPE_Z <- rpos_lognormal(n_sim, VPE_XZ_hat, SE_VPE_XZ)

draw_rg_XZ <- rcor_fisher(n_sim, rg_XZ_hat, SE_rg_XZ)
draw_rg_YZ <- rcor_fisher(n_sim, rg_YZ_hat, SE_rg_YZ)

# -----------------------------
# 8) One simulation draw -> key outputs
# -----------------------------
one_sim <- function(i){
  # --- 2-trait matrices for frontier ---
  VA_Xi <- draw_VA_X[i]; VA_Yi <- draw_VA_Y[i]
  VE_Xi <- draw_VE_X[i]; VE_Yi <- draw_VE_Y[i]
  VPE_i <- draw_VPE_XY[i]
  rg_i  <- draw_rg_XY[i]
  re_i  <- draw_re_XY[i]
  
  CovA_XY <- cov_from_r(VA_Xi, VA_Yi, rg_i)
  CovE_XY <- cov_from_r(VE_Xi, VE_Yi, re_i)
  CovP_XY <- CovA_XY + CovE_XY
  
  VP_Xi <- VA_Xi + VE_Xi + VPE_i
  VP_Yi <- VA_Yi + VE_Yi + VPE_i
  
  G2 <- matrix(c(VA_Xi, CovA_XY,
                 CovA_XY, VA_Yi),
               nrow=2, byrow=TRUE, dimnames=list(c(X,Y), c(X,Y)))
  P2 <- matrix(c(VP_Xi, CovP_XY,
                 CovP_XY, VP_Yi),
               nrow=2, byrow=TRUE, dimnames=list(c(X,Y), c(X,Y)))
  
  G2 <- make_PD(G2); P2 <- make_PD(P2)
  
  # Grid responses (CH4, MBW)
  out <- calc_resp_2trait(G2, P2, A_grid)
  ch4 <- as.numeric(out$R[1,])
  mbw <- as.numeric(out$R[2,])
  
  # Ratio point
  out_ratio <- calc_resp_2trait(G2, P2, a_ratio)
  ratio_ch4 <- as.numeric(out_ratio$R[1,1])
  ratio_mbw <- as.numeric(out_ratio$R[2,1])
  
  # Phenotypic residual point (beta_P)
  beta_P <- as.numeric(P2[X,Y] / P2[Y,Y])
  a_res <- matrix(c(-1, beta_P), ncol=1); rownames(a_res) <- c(X,Y)
  out_res <- calc_resp_2trait(G2, P2, a_res)
  resid_ch4 <- as.numeric(out_res$R[1,1])
  resid_mbw <- as.numeric(out_res$R[2,1])
  
  # ── NEW: Genetic residual point (beta_G) ──────────────────────────────────
  # beta_G = CovA(X,Y) / VA(Y) varies across MC draws because both CovA_XY
  # and VA_Y are sampled. This propagates uncertainty in the genetic
  # regression coefficient through to the response, giving wider CIs for the
  # genetic residual than the phenotypic residual if genetic parameters are
  # estimated with more uncertainty than phenotypic ones.
  beta_G <- as.numeric(G2[X,Y] / G2[Y,Y])
  a_gres <- matrix(c(-1, beta_G), ncol=1); rownames(a_gres) <- c(X,Y)
  out_gres <- calc_resp_2trait(G2, P2, a_gres)
  gresid_ch4 <- as.numeric(out_gres$R[1,1])
  gresid_mbw <- as.numeric(out_gres$R[2,1])
  
  # --- 3-trait G3 for ADG response ---
  VA_Zi <- draw_VA_Z[i]
  rg_XZ_i <- draw_rg_XZ[i]
  rg_YZ_i <- draw_rg_YZ[i]
  
  CovA_XZ <- cov_from_r(VA_Xi, VA_Zi, rg_XZ_i)
  CovA_YZ <- cov_from_r(VA_Yi, VA_Zi, rg_YZ_i)
  
  G3 <- matrix(0, 3, 3, dimnames=list(c(X,Y,Z), c(X,Y,Z)))
  G3[X,X] <- VA_Xi
  G3[Y,Y] <- VA_Yi
  G3[Z,Z] <- VA_Zi
  G3[X,Y] <- CovA_XY; G3[Y,X] <- CovA_XY
  G3[X,Z] <- CovA_XZ; G3[Z,X] <- CovA_XZ
  G3[Y,Z] <- CovA_YZ; G3[Z,Y] <- CovA_YZ
  G3 <- make_PD(G3)
  
  G_HI <- G3[, c(X,Y), drop=FALSE]
  RH <- calc_resp_H_from_B(G_HI, out$B, out$denom)
  adg <- as.numeric(RH[3,])
  
  # ADG at ratio
  RH_ratio <- calc_resp_H_from_B(G_HI, out_ratio$B, out_ratio$denom)
  ratio_adg <- as.numeric(RH_ratio[3,1])
  
  # ADG at phenotypic residual
  RH_res <- calc_resp_H_from_B(G_HI, out_res$B, out_res$denom)
  resid_adg <- as.numeric(RH_res[3,1])
  
  # ── NEW: ADG at genetic residual ──────────────────────────────────────────
  RH_gres <- calc_resp_H_from_B(G_HI, out_gres$B, out_gres$denom)
  gresid_adg <- as.numeric(RH_gres[3,1])
  
  # Lower-right quadrant summaries
  keep_lr <- (mbw > 0) & (ch4 < 0)
  if(!any(keep_lr)){
    adg_min <- adg_max <- adg_p05 <- adg_p50 <- adg_p95 <- NA_real_
    best_ch4_gok <- best_mbw_gok <- best_adg_gok <- NA_real_
  } else {
    adg_lr <- adg[keep_lr]
    adg_min <- min(adg_lr, na.rm=TRUE)
    adg_max <- max(adg_lr, na.rm=TRUE)
    adg_p05 <- as.numeric(quantile(adg_lr, 0.05, na.rm=TRUE))
    adg_p50 <- as.numeric(quantile(adg_lr, 0.50, na.rm=TRUE))
    adg_p95 <- as.numeric(quantile(adg_lr, 0.95, na.rm=TRUE))
    
    keep_gok <- keep_lr & (adg >= 0)
    if(any(keep_gok)){
      idx <- which(keep_gok)[ which.min(ch4[keep_gok]) ]
      best_ch4_gok <- ch4[idx]
      best_mbw_gok <- mbw[idx]
      best_adg_gok <- adg[idx]
    } else {
      best_ch4_gok <- best_mbw_gok <- best_adg_gok <- NA_real_
    }
  }
  
  tibble(
    ratio_ch4=ratio_ch4, ratio_mbw=ratio_mbw, ratio_adg=ratio_adg,
    resid_ch4=resid_ch4, resid_mbw=resid_mbw, resid_adg=resid_adg,
    # NEW columns for genetic residual
    gresid_ch4=gresid_ch4, gresid_mbw=gresid_mbw, gresid_adg=gresid_adg,
    adg_min=adg_min, adg_max=adg_max, adg_p05=adg_p05, adg_p50=adg_p50, adg_p95=adg_p95,
    best_ch4_gok=best_ch4_gok, best_mbw_gok=best_mbw_gok, best_adg_gok=best_adg_gok
  )
}

sim_res <- bind_rows(lapply(seq_len(n_sim), one_sim))

# -----------------------------
# 9) Results: Monte Carlo uncertainty (SD + 95% CI)
# -----------------------------
summary_table <- bind_rows(
  mc_summ(sim_res$ratio_ch4) %>% mutate(metric="ratio_ch4"),
  mc_summ(sim_res$ratio_mbw) %>% mutate(metric="ratio_mbw"),
  mc_summ(sim_res$ratio_adg) %>% mutate(metric="ratio_adg"),
  
  mc_summ(sim_res$resid_ch4) %>% mutate(metric="resid_ch4"),
  mc_summ(sim_res$resid_mbw) %>% mutate(metric="resid_mbw"),
  mc_summ(sim_res$resid_adg) %>% mutate(metric="resid_adg"),
  
  # NEW: genetic residual MC summaries
  mc_summ(sim_res$gresid_ch4) %>% mutate(metric="gresid_ch4"),
  mc_summ(sim_res$gresid_mbw) %>% mutate(metric="gresid_mbw"),
  mc_summ(sim_res$gresid_adg) %>% mutate(metric="gresid_adg"),
  
  mc_summ(sim_res$adg_min) %>% mutate(metric="adg_min_LR"),
  mc_summ(sim_res$adg_max) %>% mutate(metric="adg_max_LR"),
  mc_summ(sim_res$adg_p05) %>% mutate(metric="adg_p05_LR"),
  mc_summ(sim_res$adg_p50) %>% mutate(metric="adg_p50_LR"),
  mc_summ(sim_res$adg_p95) %>% mutate(metric="adg_p95_LR"),
  
  mc_summ(sim_res$best_ch4_gok) %>% mutate(metric="best_ch4_with_ADGge0"),
  mc_summ(sim_res$best_mbw_gok) %>% mutate(metric="best_mbw_at_bestCH4_with_ADGge0"),
  mc_summ(sim_res$best_adg_gok) %>% mutate(metric="best_adg_at_bestCH4_with_ADGge0")
) %>% select(metric, everything())

cat("\n=============================\n")
cat("MONTE CARLO SUMMARY TABLE\n")
cat("=============================\n")
print(summary_table)

# -----------------------------
# 9b) Write supplementary outputs to file
# -----------------------------
write.csv(
  summary_table,
  file = "selection_index_P3/outputs/supplementary_monte_carlo_uncertainty.csv",
  row.names = FALSE
)

full_response_grid <- resp_all_adg %>%
  select(w_CH4, w_MBW, ch4, mbw, adg) %>%
  rename(
    weight_CH4   = w_CH4,
    weight_MBW   = w_MBW,
    response_CH4 = ch4,
    response_MBW = mbw,
    response_ADG = adg
  ) %>%
  mutate(strategy = "Grid point")

# Compute ADG at the genetic residual special point for output
gresid_adg_hat <- as.numeric(calc_resp_H_from_B(G_HI_hat, out_gres_hat$B, out_gres_hat$denom)[3,1]) * 1000

special_rows <- tibble(
  weight_CH4   = c(as.numeric(a_ratio[1,1]), -1, -1),
  weight_MBW   = c(as.numeric(a_ratio[2,1]), beta_P_hat, beta_G_hat),
  response_CH4 = c(ratio_point$ch4,  resid_point$ch4,  gresid_point$ch4),
  response_MBW = c(ratio_point$mbw,  resid_point$mbw,  gresid_point$mbw),
  response_ADG = c(
    as.numeric(calc_resp_H_from_B(G_HI_hat, out_ratio_hat$B, out_ratio_hat$denom)[3,1]) * 1000,
    as.numeric(calc_resp_H_from_B(G_HI_hat, out_res_hat$B,   out_res_hat$denom)[3,1])   * 1000,
    gresid_adg_hat
  ),
  strategy = c(
    "Ratio (Taylor linearisation)",
    "Residual-P (CH4 - beta_P*MBW)",
    "Residual-G (CH4 - beta_G*MBW)"  # NEW
  )
)

full_response_grid <- bind_rows(full_response_grid, special_rows)

write.csv(
  full_response_grid,
  file = "selection_index_P3/outputs/supplementary_full_index_responses.csv",
  row.names = FALSE
)
cat("\nSupplementary files written:\n")
cat(" - supplementary_monte_carlo_uncertainty.csv\n")
cat(" - supplementary_full_index_responses.csv\n")

# -----------------------------
# 10) Error bars (95% CI) for the special points
# -----------------------------
ratio_ci <- list(
  ch4 = mc_summ(sim_res$ratio_ch4),
  mbw = mc_summ(sim_res$ratio_mbw)
)
resid_ci <- list(
  ch4 = mc_summ(sim_res$resid_ch4),
  mbw = mc_summ(sim_res$resid_mbw)
)
# NEW: genetic residual CI
gresid_ci <- list(
  ch4 = mc_summ(sim_res$gresid_ch4),
  mbw = mc_summ(sim_res$gresid_mbw)
)

pt_ci <- bind_rows(
  tibble(point="Ratio",
         x=ratio_ci$mbw$mean, y=ratio_ci$ch4$mean,
         x_l=ratio_ci$mbw$lwr_2.5, x_u=ratio_ci$mbw$upr_97.5,
         y_l=ratio_ci$ch4$lwr_2.5, y_u=ratio_ci$ch4$upr_97.5),
  tibble(point="Residual-P",
         x=resid_ci$mbw$mean, y=resid_ci$ch4$mean,
         x_l=resid_ci$mbw$lwr_2.5, x_u=resid_ci$mbw$upr_97.5,
         y_l=resid_ci$ch4$lwr_2.5, y_u=resid_ci$ch4$upr_97.5),
  # NEW row
  tibble(point="Residual-G",
         x=gresid_ci$mbw$mean, y=gresid_ci$ch4$mean,
         x_l=gresid_ci$mbw$lwr_2.5, x_u=gresid_ci$mbw$upr_97.5,
         y_l=gresid_ci$ch4$lwr_2.5, y_u=gresid_ci$ch4$upr_97.5)
)

# ============================================================
# 16) RESULTS-SUMMARY BLOCK
# ============================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

lr <- lr_hat

# Baseline special points table (point estimates, not MC)
# Now includes the genetic residual row
special_adg <- tibble(
  point = c(
    "Ratio (Taylor)",
    paste0("Residual-P (beta_P=", signif(beta_P_hat,4), ")"),
    paste0("Residual-G (beta_G=", signif(beta_G_hat,4), ")")   # NEW
  ),
  ch4_response = c(ratio_point$ch4, resid_point$ch4, gresid_point$ch4),
  mbw_response = c(ratio_point$mbw, resid_point$mbw, gresid_point$mbw),
  adg_response = c(
    as.numeric(calc_resp_H_from_B(G_HI_hat, out_ratio_hat$B, out_ratio_hat$denom)[3,1]) * 1000,
    as.numeric(calc_resp_H_from_B(G_HI_hat, out_res_hat$B,   out_res_hat$denom)[3,1])   * 1000,
    gresid_adg_hat                                                                          # NEW
  )
)

stopifnot(exists("resp_all_adg"), exists("lr"), exists("ratio_point"),
          exists("resid_point"), exists("gresid_point"), exists("special_adg"))
stopifnot(all(c("ch4","mbw","adg","w_CH4","w_MBW") %in% names(resp_all_adg)))

has_mc <- exists("summary_table") && all(c("metric","mean","mc_sd","lwr_2.5","upr_97.5","n") %in% names(summary_table))

get_mc <- function(metric_name){
  if(!has_mc) return(tibble(mean=NA_real_, mc_sd=NA_real_, lwr_2.5=NA_real_, upr_97.5=NA_real_, n=NA_integer_))
  out <- summary_table %>% filter(metric == metric_name)
  if(nrow(out) == 0) return(tibble(mean=NA_real_, mc_sd=NA_real_, lwr_2.5=NA_real_, upr_97.5=NA_real_, n=NA_integer_))
  out %>% select(mean, mc_sd, lwr_2.5, upr_97.5, n)
}

fmt_ci <- function(mu, sd, lwr, upr, digits=3){
  if(any(is.na(c(mu, sd, lwr, upr)))) return(NA_character_)
  paste0(
    signif(mu, digits), " (MC SD ", signif(sd, digits),
    "; ", signif(lwr, digits), " to ", signif(upr, digits), ")"
  )
}

# -----------------------------
# 16.1 Core points table (Ratio vs Residual-P vs Residual-G) + MC uncertainty
# -----------------------------
special_points <- special_adg %>%
  transmute(point, ch4=ch4_response, mbw=mbw_response, adg=adg_response)

cat("\n=============================\n")
cat("CORE COMPARISON: Ratio vs Residual-P vs Residual-G (point estimates)\n")
cat("=============================\n")
print(special_points)

ratio_ch4_mc  <- get_mc("ratio_ch4")
ratio_mbw_mc  <- get_mc("ratio_mbw")
ratio_adg_mc  <- get_mc("ratio_adg")
resid_ch4_mc  <- get_mc("resid_ch4")
resid_mbw_mc  <- get_mc("resid_mbw")
resid_adg_mc  <- get_mc("resid_adg")
gresid_ch4_mc <- get_mc("gresid_ch4")   # NEW
gresid_mbw_mc <- get_mc("gresid_mbw")   # NEW
gresid_adg_mc <- get_mc("gresid_adg")   # NEW

special_points_mc <- tibble(
  point = c("Ratio (Decrease CH4/MBW)", "Residual-P (X - beta_P*Y)", "Residual-G (X - beta_G*Y)"),
  ch4 = c(
    fmt_ci(ratio_ch4_mc$mean,  ratio_ch4_mc$mc_sd,  ratio_ch4_mc$lwr_2.5,  ratio_ch4_mc$upr_97.5),
    fmt_ci(resid_ch4_mc$mean,  resid_ch4_mc$mc_sd,  resid_ch4_mc$lwr_2.5,  resid_ch4_mc$upr_97.5),
    fmt_ci(gresid_ch4_mc$mean, gresid_ch4_mc$mc_sd, gresid_ch4_mc$lwr_2.5, gresid_ch4_mc$upr_97.5)
  ),
  mbw = c(
    fmt_ci(ratio_mbw_mc$mean,  ratio_mbw_mc$mc_sd,  ratio_mbw_mc$lwr_2.5,  ratio_mbw_mc$upr_97.5),
    fmt_ci(resid_mbw_mc$mean,  resid_mbw_mc$mc_sd,  resid_mbw_mc$lwr_2.5,  resid_mbw_mc$upr_97.5),
    fmt_ci(gresid_mbw_mc$mean, gresid_mbw_mc$mc_sd, gresid_mbw_mc$lwr_2.5, gresid_mbw_mc$upr_97.5)
  ),
  adg = c(
    fmt_ci(ratio_adg_mc$mean,  ratio_adg_mc$mc_sd,  ratio_adg_mc$lwr_2.5,  ratio_adg_mc$upr_97.5),
    fmt_ci(resid_adg_mc$mean,  resid_adg_mc$mc_sd,  resid_adg_mc$lwr_2.5,  resid_adg_mc$upr_97.5),
    fmt_ci(gresid_adg_mc$mean, gresid_adg_mc$mc_sd, gresid_adg_mc$lwr_2.5, gresid_adg_mc$upr_97.5)
  ),
  n = c(ratio_ch4_mc$n, resid_ch4_mc$n, gresid_ch4_mc$n)
)

cat("\n=============================\n")
cat("CORE COMPARISON: Ratio vs Residual-P vs Residual-G (MC mean + uncertainty)\n")
cat("=============================\n")
if(has_mc) {
  print(special_points_mc)
} else {
  cat("summary_table not found -> MC uncertainty not printed.\n")
}

# -----------------------------
# 16.2 Feasible ADG range in the favourable quadrant
# -----------------------------
adg_range_lr <- lr %>%
  summarise(
    n = n(),
    adg_min = min(adg, na.rm=TRUE),
    adg_max = max(adg, na.rm=TRUE),
    adg_p05 = quantile(adg, 0.05, na.rm=TRUE),
    adg_p50 = quantile(adg, 0.50, na.rm=TRUE),
    adg_p95 = quantile(adg, 0.95, na.rm=TRUE)
  )

cat("\n=============================\n")
cat("ADG RANGE UNDER CH4–MBW INDEX (favourable quadrant; point estimates)\n")
cat("=============================\n")
print(adg_range_lr)

adg_min_mc <- get_mc("adg_min_LR")
adg_max_mc <- get_mc("adg_max_LR")
adg_p05_mc <- get_mc("adg_p05_LR")
adg_p50_mc <- get_mc("adg_p50_LR")
adg_p95_mc <- get_mc("adg_p95_LR")

adg_range_lr_mc <- tibble(
  stat = c("adg_min", "adg_max", "adg_p05", "adg_p50", "adg_p95"),
  mc_summary = c(
    fmt_ci(adg_min_mc$mean, adg_min_mc$mc_sd, adg_min_mc$lwr_2.5, adg_min_mc$upr_97.5),
    fmt_ci(adg_max_mc$mean, adg_max_mc$mc_sd, adg_max_mc$lwr_2.5, adg_max_mc$upr_97.5),
    fmt_ci(adg_p05_mc$mean, adg_p05_mc$mc_sd, adg_p05_mc$lwr_2.5, adg_p05_mc$upr_97.5),
    fmt_ci(adg_p50_mc$mean, adg_p50_mc$mc_sd, adg_p50_mc$lwr_2.5, adg_p50_mc$upr_97.5),
    fmt_ci(adg_p95_mc$mean, adg_p95_mc$mc_sd, adg_p95_mc$lwr_2.5, adg_p95_mc$upr_97.5)
  ),
  n = c(adg_min_mc$n, adg_max_mc$n, adg_p05_mc$n, adg_p50_mc$n, adg_p95_mc$n)
)

cat("\n=============================\n")
cat("ADG RANGE UNDER CH4–MBW INDEX (favourable quadrant; MC mean + uncertainty)\n")
cat("=============================\n")
if(has_mc) {
  print(adg_range_lr_mc)
} else {
  cat("summary_table not found -> MC uncertainty not printed.\n")
}

# -----------------------------
# 16.3 Growth-protected region: ADG >= 0 while (ΔMBW>0, ΔCH4<0)
# -----------------------------
lr_growth_ok <- lr %>% filter(adg >= 0)

cat("\n=============================\n")
cat("GROWTH-PROTECTED REGION (ADG >= 0 within favourable quadrant)\n")
cat("=============================\n")
cat("Count:", nrow(lr_growth_ok), "out of", nrow(lr), "\n")

if(nrow(lr_growth_ok) > 0){
  best_methane_while_growth <- lr_growth_ok %>%
    arrange(ch4) %>%
    slice(1) %>%
    select(ch4, mbw, adg, w_CH4, w_MBW)
  
  cat("\nMax methane reduction achievable with ADG >= 0 (within grid):\n")
  print(best_methane_while_growth)
  
  best_growth_while_methane <- lr_growth_ok %>%
    arrange(desc(adg)) %>%
    slice(1) %>%
    select(ch4, mbw, adg, w_CH4, w_MBW)
  
  cat("\nMax ADG achievable while still in favourable quadrant and ADG >= 0:\n")
  print(best_growth_while_methane)
} else {
  cat("No growth-protected points found on this grid.\n")
}

best_ch4_mc  <- get_mc("best_ch4_with_ADGge0")
best_mbw_mc  <- get_mc("best_mbw_at_bestCH4_with_ADGge0")
best_adg_mc  <- get_mc("best_adg_at_bestCH4_with_ADGge0")

growth_protected_mc <- tibble(
  metric = c("best_ch4_with_ADGge0", "best_mbw_at_bestCH4_with_ADGge0", "best_adg_at_bestCH4_with_ADGge0"),
  mc_summary = c(
    fmt_ci(best_ch4_mc$mean, best_ch4_mc$mc_sd, best_ch4_mc$lwr_2.5, best_ch4_mc$upr_97.5),
    fmt_ci(best_mbw_mc$mean, best_mbw_mc$mc_sd, best_mbw_mc$lwr_2.5, best_mbw_mc$upr_97.5),
    fmt_ci(best_adg_mc$mean, best_adg_mc$mc_sd, best_adg_mc$lwr_2.5, best_adg_mc$upr_97.5)
  ),
  n = c(best_ch4_mc$n, best_mbw_mc$n, best_adg_mc$n)
)

cat("\n=============================\n")
cat("GROWTH-PROTECTED METRICS (MC mean + uncertainty)\n")
cat("=============================\n")
if(has_mc) {
  print(growth_protected_mc)
} else {
  cat("summary_table not found -> MC uncertainty not printed.\n")
}

# -----------------------------
# 16.5 Dominance checks (point-estimate geometry only)
# Now checks all three special points
# -----------------------------
dominates <- function(df, target_row, ch4_col="ch4", mbw_col="mbw", adg_col="adg"){
  df %>%
    filter(
      .data[[ch4_col]] <= target_row[[ch4_col]],
      .data[[mbw_col]] >= target_row[[mbw_col]],
      .data[[adg_col]] >= target_row[[adg_col]]
    ) %>%
    mutate(
      strict_better = (.data[[ch4_col]] < target_row[[ch4_col]]) |
        (.data[[mbw_col]] > target_row[[mbw_col]]) |
        (.data[[adg_col]] > target_row[[adg_col]])
    ) %>%
    filter(strict_better)
}

ratio_row  <- tibble(ch4=ratio_point$ch4,  mbw=ratio_point$mbw,  adg=special_adg$adg_response[1])
resid_row  <- tibble(ch4=resid_point$ch4,  mbw=resid_point$mbw,  adg=special_adg$adg_response[2])
gresid_row <- tibble(ch4=gresid_point$ch4, mbw=gresid_point$mbw, adg=special_adg$adg_response[3])  # NEW

dom_ratio  <- dominates(lr, ratio_row)
dom_resid  <- dominates(lr, resid_row)
dom_gresid <- dominates(lr, gresid_row)  # NEW

cat("\n=============================\n")
cat("DOMINANCE CHECKS (within favourable quadrant)\n")
cat("=============================\n")
cat("Points that dominate RATIO:", nrow(dom_ratio), "\n")
if(nrow(dom_ratio) > 0){
  cat("Example dominator of ratio:\n")
  print(dom_ratio %>% arrange(desc(adg), desc(mbw), ch4) %>% slice(1) %>%
          select(ch4, mbw, adg, w_CH4, w_MBW))
} else {
  cat("No linear-index point strictly dominates the ratio point on (CH4, MBW, ADG) simultaneously.\n")
}

cat("\nPoints that dominate RESIDUAL-P:", nrow(dom_resid), "\n")
if(nrow(dom_resid) > 0){
  cat("Example dominator of residual-P:\n")
  print(dom_resid %>% arrange(desc(adg), desc(mbw), ch4) %>% slice(1) %>%
          select(ch4, mbw, adg, w_CH4, w_MBW))
} else {
  cat("No linear-index point strictly dominates the residual-P point on (CH4, MBW, ADG) simultaneously.\n")
}

# NEW dominance check for genetic residual
cat("\nPoints that dominate RESIDUAL-G:", nrow(dom_gresid), "\n")
if(nrow(dom_gresid) > 0){
  cat("Example dominator of residual-G:\n")
  print(dom_gresid %>% arrange(desc(adg), desc(mbw), ch4) %>% slice(1) %>%
          select(ch4, mbw, adg, w_CH4, w_MBW))
} else {
  cat("No linear-index point strictly dominates the residual-G point on (CH4, MBW, ADG) simultaneously.\n")
}

# -----------------------------
# 16.6 Index flexibility advantage
# -----------------------------
cat("\n=============================\n")
cat("INDEX FLEXIBILITY: ADG span vs fixed points\n")
cat("=============================\n")
cat("ADG range in favourable quadrant: [",
    signif(adg_range_lr$adg_min, 4), ", ", signif(adg_range_lr$adg_max, 4), "]\n", sep="")
cat("ADG at ratio:      ", signif(special_adg$adg_response[1], 4), "\n", sep="")
cat("ADG at residual-P: ", signif(special_adg$adg_response[2], 4), "\n", sep="")
cat("ADG at residual-G: ", signif(special_adg$adg_response[3], 4), "\n", sep="")

# -----------------------------
# 16.8 Compact summary table for Results (point estimates)
# -----------------------------
results_compact <- bind_rows(
  special_points %>% mutate(section="Special points"),
  tibble(point="Linear index (favourable quadrant)",
         ch4=NA_real_, mbw=NA_real_, adg=NA_real_, section="Feasible region") %>%
    mutate(
      ch4 = min(lr$ch4, na.rm=TRUE),
      mbw = max(lr$mbw, na.rm=TRUE),
      adg = max(lr$adg, na.rm=TRUE)
    )
) %>%
  select(section, point, ch4, mbw, adg)

cat("\n=============================\n")
cat("COMPACT SUMMARY TABLE (point estimates)\n")
cat("=============================\n")
print(results_compact)

if(has_mc){
  mc_compact <- tibble(
    metric = c("ratio_ch4","ratio_mbw","ratio_adg",
               "resid_ch4","resid_mbw","resid_adg",
               "gresid_ch4","gresid_mbw","gresid_adg",   # NEW
               "adg_min_LR","adg_p05_LR","adg_p50_LR","adg_p95_LR","adg_max_LR",
               "best_ch4_with_ADGge0","best_mbw_at_bestCH4_with_ADGge0","best_adg_at_bestCH4_with_ADGge0")
  ) %>%
    rowwise() %>%
    mutate(
      mean     = get_mc(metric)$mean,
      mc_sd    = get_mc(metric)$mc_sd,
      lwr_2.5  = get_mc(metric)$lwr_2.5,
      upr_97.5 = get_mc(metric)$upr_97.5,
      n        = get_mc(metric)$n
    ) %>% ungroup()
  
  cat("\n=============================\n")
  cat("MONTE CARLO SUMMARY TABLE (for reporting)\n")
  cat("=============================\n")
  print(mc_compact)
} else {
  cat("\n(MC summary_table not found -> skipping MC compact table.)\n")
}

# ============================================================
# FRONTIER PLOT v5 — HULL LINE COLOURED BY ADG RESPONSE
# Now shows THREE special points: ratio, residual-P, residual-G
# ============================================================
suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ggtext)
})

lr <- resp_all_adg %>% filter(mbw > 0, ch4 < 0)

hull_lr_adg <- hull_df %>%
  filter(mbw > 0, ch4 < 0) %>%
  mutate(mbw_r = round(mbw, 6), ch4_r = round(ch4, 6)) %>%
  left_join(
    resp_all_adg %>%
      mutate(mbw_r = round(mbw, 6), ch4_r = round(ch4, 6)) %>%
      select(mbw_r, ch4_r, adg),
    by = c("mbw_r", "ch4_r")
  )

if (any(is.na(hull_lr_adg$adg))) {
  missing_idx <- which(is.na(hull_lr_adg$adg))
  for (i in missing_idx) {
    dists <- sqrt(
      (resp_all_adg$mbw - hull_lr_adg$mbw[i])^2 +
        (resp_all_adg$ch4 - hull_lr_adg$ch4[i])^2
    )
    hull_lr_adg$adg[i] <- resp_all_adg$adg[which.min(dists)]
  }
}

adg_lims        <- range(lr$adg, na.rm = TRUE)
adg_pad         <- diff(adg_lims) * 0.05
adg_lims_padded <- c(adg_lims[1] - adg_pad, adg_lims[2] + adg_pad)
pct_positive    <- round(100 * mean(lr$adg >= 0, na.rm = TRUE), 1)
med_adg         <- median(lr$adg, na.rm = TRUE)

palette_stops <- c("#8B0000", "#CC3300", "#FF6644", "#FAEBD7",
                   "#88CC88", "#2E8B2E", "#145214")

value_stops <- rescale(c(
  adg_lims[1], adg_lims[1]*0.5, adg_lims[1]*0.1, 0,
  adg_lims[2]*0.1, adg_lims[2]*0.5, adg_lims[2]
))

adg_col_scale <- scale_colour_gradientn(
  colours = palette_stops, values = value_stops, limits = adg_lims,
  breaks  = c(adg_lims[1], adg_lims[2]),
  labels  = sprintf("%.2f", c(adg_lims[1], adg_lims[2])),
  name    = expression("ADG correlated response (g day"^{-1}*")")
)

p_frontier <- ggplot() +
  
  # Background cloud
  geom_point(
    data  = lr,
    aes(x = mbw, y = ch4, colour = adg),
    alpha = 0.25, size = 1.4
  ) +
  adg_col_scale +
  
  # Hull line coloured by ADG
  geom_path(
    data      = hull_lr_adg,
    aes(x = mbw, y = ch4, colour = adg),
    linewidth = 2.2, lineend = "round", linejoin = "round"
  ) +
  
  # Ratio point
  geom_point(
    data = ratio_point, aes(x = mbw, y = ch4),
    size = 5, shape = 24, fill = "white", colour = "#1A1A2E", stroke = 1.5
  ) +
  annotate(
    "richtext",
    x     = ratio_point$mbw + 0.02, y = ratio_point$ch4 - 0.05,
    label = paste0(
      "<b>Ratio</b><br>",
      "\u0394CH4: ", round(ratio_point$ch4, 2), " (\u00b1", round(get_mc("ratio_ch4")$mc_sd, 2), ")<br>",
      "\u0394MBW: ", round(ratio_point$mbw, 2),  " (\u00b1", round(get_mc("ratio_mbw")$mc_sd, 2), ")<br>",
      "\u0394ADG: ", round(special_adg$adg_response[1], 2), " (\u00b1", round(get_mc("ratio_adg")$mc_sd * 1000, 2), ")"
    ),
    size = 3.3, colour = "#1A1A2E", hjust = 0, fill = NA, label.color = NA
  ) +
  
  # Phenotypic residual point (square, unchanged)
  geom_point(
    data = resid_point %>% filter(mbw > 0, ch4 < 0),
    aes(x = mbw, y = ch4),
    size = 5, shape = 22, fill = "white", colour = "#1A1A2E", stroke = 1.5
  ) +
  annotate(
    "richtext",
    x     = resid_point$mbw + 0.02, y = resid_point$ch4 - 0.05,
    label = paste0(
      "<b>Residual-P</b><br>",
      "\u0394CH4: ", round(resid_point$ch4, 2), " (\u00b1", round(get_mc("resid_ch4")$mc_sd, 2), ")<br>",
      "\u0394MBW: ", round(resid_point$mbw, 2),  " (\u00b1", round(get_mc("resid_mbw")$mc_sd, 2), ")<br>",
      "\u0394ADG: ", round(special_adg$adg_response[2], 2), " (\u00b1", round(get_mc("resid_adg")$mc_sd * 1000, 2), ")"
    ),
    size = 3.3, colour = "#1A1A2E", hjust = 0, fill = NA, label.color = NA
  ) +
  
  # NEW: Genetic residual point (diamond)
  geom_point(
    data = gresid_point %>% filter(mbw > 0, ch4 < 0),
    aes(x = mbw, y = ch4),
    size = 5, shape = 23, fill = "white", colour = "#1A1A2E", stroke = 1.5
  ) +
  annotate(
    "richtext",
    x     = gresid_point$mbw + 0.02, y = gresid_point$ch4 - 0.05,
    label = paste0(
      "<b>Residual-G</b><br>",
      "\u0394CH4: ", round(gresid_point$ch4, 2), " (\u00b1", round(get_mc("gresid_ch4")$mc_sd, 2), ")<br>",
      "\u0394MBW: ", round(gresid_point$mbw, 2),  " (\u00b1", round(get_mc("gresid_mbw")$mc_sd, 2), ")<br>",
      "\u0394ADG: ", round(special_adg$adg_response[3], 2), " (\u00b1", round(get_mc("gresid_adg")$mc_sd * 1000, 2), ")"
    ),
    size = 3.3, colour = "#1A1A2E", hjust = 0, fill = NA, label.color = NA
  ) +
  
  # Reference lines
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.35) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.35) +
  
  labs(
    x = expression("Predicted response in metabolic body weight (kg)"),
    y = expression("Predicted response in CH"[4]*" production (g day"^{-1}*")")
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_blank(),
    plot.subtitle    = element_blank(),
    axis.title       = element_text(size = 10, colour = "#1A1A2E"),
    legend.position  = "bottom",
    legend.title     = element_text(size = 9),
    legend.key.width = unit(2.5, "cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey93")
  )

p_frontier

ggsave(
  "selection_index_P3/outputs/frontier_adg_v5_with_genetic_residual.png",
  width  = 8,
  height = 8,
  dpi    = 300,
  bg     = "white"
)