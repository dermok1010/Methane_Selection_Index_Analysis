
# ============================================================
# SELECTION RESPONSE COMPARISON OF METHANE TRAITS WITH MONTE CARLO UNCERTAINTY
# Uses ASReml .pvc estimates (estimate followed by SE)
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
VA_Y_hat  <- 1.51857   ; SE_VA_Y  <- 0.164883    # Additive genetic variance: MBW (anchors G2 and G3 diagonal)
VA_X_hat  <- 2.72527   ; SE_VA_X  <- 0.251641    # Additive genetic variance: CH4 (anchors G2 and G3 diagonal)
rg_XY_hat <- -0.033    ; SE_rg_XY <- 0.0838      # Genetic correlation: MBW-CH4 (near-zero; weak genetic association)

# Residual
VE_Y_hat  <- 0.964220  ; SE_VE_Y  <- 0.168305E-01  # Residual variance: MBW
VE_X_hat  <- 10.7364   ; SE_VE_X  <- 0.170771      # Residual variance: CH4
re_XY_hat <- 0.2358    ; SE_re_XY <- 0.0115         # Residual correlation: MBW-CH4 (modest positive; larger than rg)

# PE: ide(ANI_ID) shared across both traits
VPE_XY_hat <- 2.38080  ; SE_VPE_XY <- 0.157983     # Permanent environment variance: ide(ANI_ID), shared scalar across MBW and CH4
# Inflates VP diagonal only; no PE covariance between traits

# ---- Bivariate: Metabolic_BW (trait1) vs adg (trait2) ----
VA_Y2_hat  <- 4.29765        ; SE_VA_Y2  <- 0.890151E-01   # Additive genetic variance: MBW (this bivariate's estimate)
VA_Z_from_Y_hat <- 0.194251E-02 ; SE_VA_Z_from_Y <- 0.194836E-03  # Additive genetic variance: ADG (from MBW-ADG bivariate)
rg_YZ_hat  <- 0.6878         ; SE_rg_YZ  <- 0.0341         # Genetic correlation: MBW-ADG
VE_Y2_hat  <- 0.994608       ; SE_VE_Y2  <- 0.174554E-01   # Residual variance: MBW (this bivariate's estimate)
VE_Z_from_Y_hat <- 0.134158E-03 ; SE_VE_Z_from_Y <- 0.516191E-05  # Residual variance: ADG (from MBW-ADG bivariate)
re_YZ_hat  <- 0.0830         ; SE_re_YZ  <- 0.0343         # Residual correlation: MBW-ADG
VPE_YZ_hat <- 0.793797E-03   ; SE_VPE_YZ <- 0.167115E-03  # Permanent environment variance: ide(ANI_ID), MBW-ADG bivariate

# ---- Bivariate: ch4_g_day2_1v3 (trait1) vs adg (trait2) ----
VA_X2_hat  <- 4.62718        ; SE_VA_X2  <- 0.230782       # Additive genetic variance: CH4 (this bivariate's estimate)
VA_Z_from_X_hat <- 0.149254E-02 ; SE_VA_Z_from_X <- 0.205584E-03  # Additive genetic variance: ADG (from CH4-ADG bivariate; used in G3 assembly)
rg_XZ_hat  <- 0.4104         ; SE_rg_XZ  <- 0.0486         # Genetic correlation: CH4-ADG (used in G3 assembly)
VE_X2_hat  <- 10.9767        ; SE_VE_X2  <- 0.176163       # Residual variance: CH4 (this bivariate's estimate)
VE_Z_from_X_hat <- 0.133262E-03 ; SE_VE_Z_from_X <- 0.510974E-05  # Residual variance: ADG (from CH4-ADG bivariate)
re_XZ_hat  <- 0.0846         ; SE_re_XZ  <- 0.0300         # Residual correlation: CH4-ADG
VPE_XZ_hat <- 0.817222E-03   ; SE_VPE_XZ <- 0.183645E-03  # Permanent environment variance: ide(ANI_ID), CH4-ADG bivariate (used in G3 assembly)



# -----------------------------
# 2) Helpers
# -----------------------------
# Returns the smallest eigenvalue of a matrix.
# Eigenvalues describe the "stretching factors" of a matrix along its principal
# axes. For a variance-covariance matrix to be valid (positive definite),
# ALL eigenvalues must be positive. A zero or negative eigenvalue means the
# ellipsoid has collapsed into a flat plane or flipped inside out which is
# physically impossible for a real covariance structure.
eig_min <- function(M) min(eigen(M, symmetric=TRUE, only.values=TRUE)$values)

# Checks if a matrix is valid (positive definite) and repairs it if not.
# This can happen in Monte Carlo draws when sampled parameters produce a
# mathematically impossible covariance structure — e.g. a genetic correlation
# that implies more shared variance than actually exists in either trait.
# nearPD() finds the nearest valid positive definite matrix by nudging the
# eigenvalues just above zero. Think of it as inflating a slightly collapsed
# ellipsoid back into a proper 3D shape.
# The threshold 1e-10 rather than 0 guards against floating point rounding
# errors producing tiny negative eigenvalues in otherwise valid matrices.
make_PD <- function(M){
  if(eig_min(M) < 1e-10) as.matrix(nearPD(M)$mat) else M
}

# ── PARAMETER SAMPLING FOR MONTE CARLO ───────────────────────────────────────

# Draws n random values for a variance component (VA, VE, VPE etc.)
# centred on the point estimate x with uncertainty governed by se.
#
# Why log-normal and not normal?
# Variance components must be strictly positive — you cannot have negative
# genetic variance. A normal distribution centred on e.g. VA=3.0 would
# occasionally produce draws of -0.5, which is nonsensical. The log-normal
# is always positive and its shape naturally reflects the right-skew of
# variance component sampling distributions.
#
# How the log-normal is parameterised to match x and se:
# cv = se/x  is the coefficient of variation — relative uncertainty
# sigma and mu are the log-scale parameters back-calculated so that the
# resulting distribution has mean = x and SD = se on the original scale.
# The exp() at the end transforms back from log-scale to the original scale.
rpos_lognormal <- function(n, x, se){
  if(is.na(se) || se <= 0) return(rep(x, n))
  cv <- se / x
  sigma <- sqrt(log(1 + cv^2))
  mu <- log(x) - 0.5 * sigma^2
  exp(rnorm(n, mean=mu, sd=sigma))
}

# Draws n random values for a correlation (rg or re)
# centred on the point estimate r with uncertainty governed by se_r.
#
# Why not just sample r directly from a normal distribution?
# Correlations are bounded between -1 and +1. A normal distribution
# centred on e.g. rg=0.07 with SE=0.077 would sometimes produce draws
# outside this range, which are impossible. Near the boundaries the
# sampling distribution of r is also skewed, not symmetric.
#
# The Fisher-z transformation solves both problems:
# z = atanh(r) transforms r from (-1,+1) to (-inf,+inf) where normal
# sampling is appropriate. The SE is also transformed to the z scale.
# After drawing on the z scale we back-transform via tanh() which
# automatically keeps results within (-1,+1).
# The pmax/pmin clipping to ±0.999 prevents exact ±1 which would cause
# downstream problems (singular matrices).
rcor_fisher <- function(n, r, se_r){
  if(is.na(se_r) || se_r <= 0) return(rep(r, n))
  r <- max(min(r, 0.999), -0.999)
  z <- atanh(r)
  se_z <- se_r / (1 - r^2)
  z_draw <- rnorm(n, mean=z, sd=se_z)
  pmax(pmin(tanh(z_draw), 0.999), -0.999)
}

# ── COVARIANCE RECONSTRUCTION ────────────────────────────────────────────────

# Reconstructs the covariance between two traits from their variances and
# their correlation. ASReml reports rg and the individual variances separately
# rather than reporting CovA directly. Each Monte Carlo draw recalculates
# the covariance from freshly sampled variances and correlation, so the
# uncertainty in the covariance inherits from all three sources simultaneously.
# This is the off-diagonal entry that goes into G or the residual part of P.
cov_from_r <- function(Vi, Vj, r){
  r * sqrt(Vi * Vj)
}
# ── INDEX WEIGHT COMPUTATION ─────────────────────────────────────────────────

# Computes the optimal index weights b = P^{-1} G a for a single desired
# direction a. solve(P, x) computes P^{-1} x without explicitly inverting P
# which is numerically more stable. G %*% a maps the desired genetic direction
# through the genetic covariance structure first, then P^{-1} rescales by the
# phenotypic structure to give the actual measurement weights to apply to
# observed animal phenotypes.
compute_b <- function(G, P, a) solve(P, G %*% a)

# ── VECTORISED RESPONSE CALCULATION ACROSS GRID ──────────────────────────────

# Computes the selection index response for ALL index weight combinations
# in the grid simultaneously.
# A is a 2 x m matrix where each column is one desired direction vector,
# and m is the total number of grid points.
#
# Step by step:
# G2 %*% A  : maps all m desired directions through genetic structure (2 x m)
# solve(P2, ...): applies P^{-1} to each column — gives B, the optimal index
#                 weights for each of the m grid directions (2 x m)
# P2 %*% B  : applies P to B — gives the phenotypic variance of each index
# colSums(B * (P2 %*% B)): computes b'Pb for each column — variance of index scores
# sqrt(...)  : gives SD of index scores = the normalising denominator (length m)
# t(G2) %*% B : applies genetic structure to index weights — raw genetic response
# dividing by denom: normalises to i=1 selection intensity
#
# Returns B (index weights), denom (normalising SDs), and R (responses)
# as a list so they can be reused without recalculation for the ADG step.
calc_resp_2trait <- function(G2, P2, A){
  # A is 2 x m
  B <- solve(P2, G2 %*% A)               # optimal index weights for all grid points
  denom <- sqrt(colSums(B * (P2 %*% B)))  # SD of index scores for each grid point
  R <- (t(G2) %*% B) / matrix(denom, nrow=2, ncol=length(denom), byrow=TRUE)
  list(B=B, denom=denom, R=R)
}

# ── CORRELATED RESPONSE IN A THIRD TRAIT ─────────────────────────────────────

# Computes the correlated genetic response in traits NOT in the index
# (here ADG) given the already-computed index weights B and normalising
# denominators from the 2-trait step above.
#
# G_HI is the submatrix of the 3-trait G containing the genetic covariances
# between ALL traits (rows) and only the INDEXED traits (columns).
# Multiplying by B gives the raw genetic response in all traits including ADG.
# Dividing by denom normalises to i=1 as before.
# ADG comes along for the ride purely through its genetic correlations with
# CH4 and MBW — it is not in the index, it is a correlated passenger.
calc_resp_H_from_B <- function(G_HI, B, denom){
  (G_HI %*% B) / matrix(denom, nrow=nrow(G_HI), ncol=length(denom), byrow=TRUE)
}

# ── FRONTIER IDENTIFICATION ───────────────────────────────────────────────────

# Finds the convex hull of the response cloud — the outer boundary of all
# achievable genetic responses across the grid of index weights.
# chull() returns the indices of points forming the convex hull in
# counter-clockwise order. These are the Pareto-optimal points — responses
# where you cannot improve one trait without sacrificing the other.
# Points inside the hull are achievable but suboptimal.
# The nrow < 3 guard handles degenerate cases where too few valid points exist.
get_hull <- function(df, x="mbw", y="ch4"){
  if(nrow(df) < 3) return(df[0, , drop=FALSE])
  idx <- chull(df[[x]], df[[y]])
  df[idx, , drop=FALSE]
}
# ── MONTE CARLO SUMMARY STATISTICS ───────────────────────────────────────────

# Summarises the distribution of a quantity across the 1000 Monte Carlo draws.
# Filters to finite values first to handle any draws that produced NA or Inf
# (e.g. from a degenerate matrix that nearPD could not fully repair).
# Returns the mean, MC standard deviation, and 95% credible interval
# (2.5th and 97.5th percentiles) — these are the uncertainty bands that
# appear in the summary table and could be plotted around the frontier.
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
# Define the range of weights to sweep for each trait.
# These are the "desired direction" vectors a — every combination of
# w_CH4 and w_MBW represents a different breeding objective.
# e.g. w_CH4 = -1, w_MBW = 0 means "only reduce methane, ignore MBW"
#      w_CH4 = -1, w_MBW = 1 means "reduce methane AND increase MBW equally"
#      w_CH4 = -5, w_MBW = 1 means "strongly prioritise methane reduction"
# The range -10 to +10 is arbitrary but wide enough to capture all
# meaningful directions. The actual scale doesn't matter because responses
# are normalised to i=1 — only the RATIO of weights affects the direction.
w_CH4_seq <- seq(-10, 10, by=0.5)
w_MBW_seq <- seq(-10, 10, by=0.5)

# Create every combination of w_CH4 and w_MBW — this is the full grid.
# expand.grid gives all pairwise combinations: 41 x 41 = 1,681 directions.
# Remove the (0,0) combination because a zero weight vector has no direction
# and would cause a division by zero in the response calculation.
grid <- expand.grid(w_CH4 = w_CH4_seq, w_MBW = w_MBW_seq) %>%
  filter(!(w_CH4 == 0 & w_MBW == 0))

# Stack the weight vectors as columns of a 2 x m matrix A.
# Each column is one desired direction vector a = [w_CH4, w_MBW].
# This vectorised form allows calc_resp_2trait to compute all 1,680
# responses simultaneously as matrix operations rather than looping.
# Rows are named after the traits so the matrix algebra is self-documenting.
A_grid <- rbind(grid$w_CH4, grid$w_MBW)
rownames(A_grid) <- c(X, Y)

# ── SPECIAL POINT 1: RATIO (TAYLOR LINEARISATION) ────────────────────────────

# The ratio objective targets methane efficiency: minimise CH4/MBW.
# Directly minimising a ratio is nonlinear so we use a first-order
# Taylor expansion around the trait means to get approximate linear weights.
# For a ratio f = X/Y, the gradient (economic weights) at the mean are:
#   df/dX =  1/mu_Y        (reducing X by 1 unit reduces ratio by 1/mu_Y)
#   df/dY = -mu_X/mu_Y^2   (increasing Y by 1 unit reduces ratio by mu_X/mu_Y^2)
# We want to MINIMISE the ratio so we negate both:
#   a_CH4 = -1/mu_Y        (penalise CH4)
#   a_MBW = +mu_X/mu_Y^2  (reward MBW, because more MBW dilutes methane)
# This is a single column matrix (one specific desired direction).
a_ratio <- matrix(c(-1/mu_Y, mu_X/(mu_Y^2)), ncol=1)
rownames(a_ratio) <- c(X, Y)

# -----------------------------
# 4) Baseline matrices (use the CH4–MBW bivariate as baseline for 2-trait frontier)
# -----------------------------
# Reconstruct the genetic covariance between CH4 and MBW from the
# point estimates of the individual variances and genetic correlation.
# This is the off-diagonal entry shared by both G and the genetic part of P.
CovA_XY_hat <- cov_from_r(VA_X_hat, VA_Y_hat, rg_XY_hat)
# Same for the residual covariance — the environmental component of
# the phenotypic covariance between traits.
CovE_XY_hat <- cov_from_r(VE_X_hat, VE_Y_hat, re_XY_hat)
# Phenotypic covariance = genetic covariance + residual covariance.
# Note: PE (permanent environment) contributes nothing to the off-diagonal
CovP_XY_hat <- CovA_XY_hat + CovE_XY_hat    # no PE covariance
# Total phenotypic variance for each trait = genetic + residual + PE.
# VPE is the same value for both traits because it is a shared scalar
# from ide(ANI_ID). It makes both ellipsoid axes larger but does not
# tilt the ellipsoid (no covariance contribution).
VP_X_hat <- VA_X_hat + VE_X_hat + VPE_XY_hat
VP_Y_hat <- VA_Y_hat + VE_Y_hat + VPE_XY_hat

# Assemble the 2x2 additive genetic variance-covariance matrix G.
# Diagonals: additive genetic variances for CH4 and MBW.
# Off-diagonals: additive genetic covariance (symmetric).
# G describes the shape of the hidden genetic ellipsoid.
# dimnames labels rows and columns with trait names for clarity.
G2_hat <- matrix(c(VA_X_hat, CovA_XY_hat,
                   CovA_XY_hat, VA_Y_hat),
                 nrow=2, byrow=TRUE, dimnames=list(c(X,Y), c(X,Y)))

# Assemble the 2x2 phenotypic variance-covariance matrix P.
# Diagonals: total phenotypic variances (VA + VE + VPE) for each trait.
# Off-diagonals: phenotypic covariance (CovA + CovE, no PE contribution).
# P describes the shape of the observable phenotypic ellipsoid —
# always larger and potentially differently tilted than G.
P2_hat <- matrix(c(VP_X_hat, CovP_XY_hat,
                   CovP_XY_hat, VP_Y_hat),
                 nrow=2, byrow=TRUE, dimnames=list(c(X,Y), c(X,Y)))

# Check both matrices are positive definite (valid ellipsoids).
# With point estimates this should always pass but good practice.
# If any eigenvalue < 1e-10 the matrix is repaired via nearPD().
G2_hat <- make_PD(G2_hat)
P2_hat <- make_PD(P2_hat)

# Run the response calculation across all 1,680 grid directions at once.
# Returns a list containing:
#   $B     — optimal index weights for each grid direction (2 x 1680 matrix)
#   $denom — normalising SD of index scores for each direction (length 1680)
#   $R     — normalised genetic responses in CH4 and MBW (2 x 1680 matrix)
out_hat <- calc_resp_2trait(G2_hat, P2_hat, A_grid)
# Attach the CH4 and MBW responses back to the grid dataframe.
# Row 1 of R is CH4 response, row 2 is MBW response for each grid point.
# Each row is now one point in the 2D genetic response space.
resp_all <- grid %>%
  mutate(ch4 = as.numeric(out_hat$R[1,]),
         mbw = as.numeric(out_hat$R[2,]))
# Find the convex hull of the response cloud — the frontier.
# These are the Pareto-optimal points: the outer boundary of all
# achievable (CH4, MBW) response combinations.
# Closing the hull by appending the first point creates a closed polygon
# for plotting with geom_path().
hull_df <- get_hull(resp_all, x="mbw", y="ch4") %>% mutate(type="Frontier")
hull_df_closed <- bind_rows(hull_df, hull_df[1, , drop=FALSE])

# Compute the response at the Taylor ratio objective specifically.
# This is a single column so R is a 2x1 matrix — one point on the frontier.
# This tells you: if you select to minimise CH4/MBW using these index weights,
# where does that land in the (CH4 response, MBW response) space?
out_ratio_hat <- calc_resp_2trait(G2_hat, P2_hat, a_ratio)
ratio_point <- tibble(type="Ratio (Taylor)",
                      ch4 = as.numeric(out_ratio_hat$R[1,1]),
                      mbw = as.numeric(out_ratio_hat$R[2,1]))

# The residual index selects on CH4 after removing the phenotypic
# regression of CH4 on MBW. beta_P is that phenotypic regression coefficient
# — how much CH4 changes per unit of MBW in the phenotypic data.
# This comes entirely from P (not G), reflecting observable associations
# rather than genetic ones.
# The index a = [-1, beta_P] says: penalise CH4 but give credit for MBW
# in proportion to how much MBW phenotypically predicts CH4.
# This targets animals that produce less methane than expected for their
# body size — size-independent methane efficiency.
beta_P_hat <- as.numeric(P2_hat[X,Y] / P2_hat[Y,Y])
a_res_hat <- matrix(c(-1, beta_P_hat), ncol=1); rownames(a_res_hat) <- c(X,Y)
out_res_hat <- calc_resp_2trait(G2_hat, P2_hat, a_res_hat)
resid_point <- tibble(type=paste0("Residual (beta_P=", signif(beta_P_hat,4), ")"),
                      ch4 = as.numeric(out_res_hat$R[1,1]),
                      mbw = as.numeric(out_res_hat$R[2,1]))

# ── SECTION 5: BUILD THE 3x3 G MATRIX FOR ADG ────────────────────────────────

# We have three pairwise bivariate models but no single 3-trait joint fit.
# We therefore assemble G3 from the pairwise estimates, choosing ONE
# source for each trait's own variance to keep the diagonal consistent:
#   CH4 (X): variance from the CH4-MBW bivariate (anchors the frontier)
#   MBW (Y): variance from the CH4-MBW bivariate (same model, consistent)
#   ADG (Z): variance from the CH4-ADG bivariate (arbitrary but documented)
# This is an approximation — a true 3-trait fit would be preferable but
# is not always computationally feasible or available from ASReml output.

VA_X3_hat <- VA_X_hat
VE_X3_hat <- VE_X_hat
VPE_X3_hat <- VPE_XY_hat

VA_Y3_hat <- VA_Y_hat
VE_Y3_hat <- VE_Y_hat
VPE_Y3_hat <- VPE_XY_hat

VA_Z3_hat <- VA_Z_from_X_hat
VE_Z3_hat <- VE_Z_from_X_hat
VPE_Z3_hat <- VPE_XZ_hat

# Reconstruct the three pairwise additive genetic covariances.
# CovA_XY already computed in section 4.
# CovA_XZ uses rg from the CH4-ADG bivariate.
# CovA_YZ uses rg from the MBW-ADG bivariate.
# Note: because these rg values come from different models, the resulting
# G3 matrix may not be perfectly positive definite — make_PD() handles this.
CovA_XZ_hat <- cov_from_r(VA_X3_hat, VA_Z3_hat, rg_XZ_hat)
CovA_YZ_hat <- cov_from_r(VA_Y3_hat, VA_Z3_hat, rg_YZ_hat)
# Build the 3x3 G matrix by slotting variances onto the diagonal
# and covariances onto the off-diagonals (symmetrically).
# Starting from a zero matrix and filling by name is safer than
# byrow= construction for a 3x3 as it makes the structure explicit.
G3_hat <- matrix(0, 3, 3, dimnames=list(c(X,Y,Z), c(X,Y,Z)))
G3_hat[X,X] <- VA_X3_hat
G3_hat[Y,Y] <- VA_Y3_hat
G3_hat[Z,Z] <- VA_Z3_hat
G3_hat[X,Y] <- CovA_XY_hat; G3_hat[Y,X] <- CovA_XY_hat
G3_hat[X,Z] <- CovA_XZ_hat; G3_hat[Z,X] <- CovA_XZ_hat
G3_hat[Y,Z] <- CovA_YZ_hat; G3_hat[Z,Y] <- CovA_YZ_hat
# Repair G3 if pairwise assembly has produced any invalid eigenvalues.
# This is more likely here than for G2 because pairwise-assembled matrices
# are not guaranteed to be globally consistent.
G3_hat <- make_PD(G3_hat)

# ── COMPUTE ADG CORRELATED RESPONSE ACROSS THE GRID ──────────────────────────

# Extract only the columns of G3 corresponding to the INDEXED traits (CH4, MBW).
# This 3x2 submatrix is what connects the 2-trait index to the 3-trait response.
# Rows = all three traits (CH4, MBW, ADG).
# Columns = only the two traits in the index.
# ADG is NOT in the index — it responds only through its genetic correlations
# with CH4 and MBW encoded in the off-diagonal entries of this submatrix.
G_HI_hat <- G3_hat[, c(X,Y), drop=FALSE]
# Compute the correlated response in all three traits across the full grid.
# Reuses the index weights B and denominators denom already computed in
# out_hat — no need to redo the P inversion.
# RH_hat is a 3 x 1680 matrix: rows are CH4, MBW, ADG responses,
# columns are the 1680 grid directions.
RH_hat <- calc_resp_H_from_B(G_HI_hat, out_hat$B, out_hat$denom)
# Attach the ADG correlated response to the existing response dataframe.
# Row 3 of RH_hat is the ADG response for each grid point.
# Each row now has a CH4 response, MBW response, and ADG correlated response —
# a complete 3D picture of what happens under each index weight combination.
resp_all_adg <- resp_all %>% mutate(adg = as.numeric(RH_hat[3,]) * 1000)
# Filter to the desirable quadrant: index directions that simultaneously
# increase MBW (mbw > 0) AND decrease CH4 (ch4 < 0).
# This subset is the region of practical interest — the part of the frontier
# where both primary breeding goals are moving in the right direction.
# Within this region, adg values show whether growth is helped or harmed.
lr_hat <- resp_all_adg %>% filter(mbw > 0, ch4 < 0)

# -----------------------------
# 6) Baseline plots (unchanged visuals)
# -----------------------------
ggplot() +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  geom_point(data=resp_all, aes(x=mbw, y=ch4), alpha=0.12, size=1) +
  geom_path(data=hull_df_closed, aes(x=mbw, y=ch4), linewidth=1.0) +
  geom_point(data=ratio_point, aes(x=mbw, y=ch4), size=4, shape=17) +
  geom_point(data=resid_point, aes(x=mbw, y=ch4), size=4, shape=15) +
  labs(
    x="Predicted response in Metabolic_BW (i=1 normalised)",
    y="Predicted response in CH4/day (i=1 normalised)",
    title="2D response region and frontier for CH4 & Metabolic_BW",
    subtitle="Line: convex-hull frontier. Triangle: ratio (Taylor). Square: residual (X−βY)."
  ) +
  theme_minimal()

ggplot() +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  geom_point(data=resp_all, aes(x=mbw, y=ch4), alpha=0.12, size=1) +
  geom_path(data=hull_df_closed, aes(x=mbw, y=ch4), linewidth=1.0) +
  geom_point(data=ratio_point, aes(x=mbw, y=ch4), size=4, shape=17) +
  geom_point(data=resid_point, aes(x=mbw, y=ch4), size=4, shape=15) +
  coord_cartesian(
    xlim = c(0, max(resp_all$mbw, na.rm=TRUE)),
    ylim = c(min(resp_all$ch4, na.rm=TRUE), 0)
  ) +
  labs(
    x="Predicted response in Metabolic_BW (i=1 normalised)",
    y="Predicted response in CH4/day (i=1 normalised)",
    title="2D response region and frontier for CH4 & Metabolic_BW",
    subtitle="Zoomed to bottom-right quadrant (↑MBW, ↓CH4)"
  ) +
  theme_minimal()

# -----------------------------
# 7) Monte Carlo
# -----------------------------
set.seed(1)
n_sim <- 1000

# Draw CH4–MBW parameters (for 2-trait frontier + beta_P)
draw_VA_X  <- rpos_lognormal(n_sim, VA_X_hat,  SE_VA_X)
draw_VA_Y  <- rpos_lognormal(n_sim, VA_Y_hat,  SE_VA_Y)
draw_VE_X  <- rpos_lognormal(n_sim, VE_X_hat,  SE_VE_X)
draw_VE_Y  <- rpos_lognormal(n_sim, VE_Y_hat,  SE_VE_Y)
draw_VPE_XY <- rpos_lognormal(n_sim, VPE_XY_hat, SE_VPE_XY)

draw_rg_XY <- rcor_fisher(n_sim, rg_XY_hat, SE_rg_XY)
draw_re_XY <- rcor_fisher(n_sim, re_XY_hat, SE_re_XY)

# Draw CH4–ADG and MBW–ADG pieces for ADG response uncertainty in G3
# ADG variances: choose from CH4–ADG bivariate for consistency with earlier assembly
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
  
  # Residual point (beta_P)
  beta_P <- as.numeric(P2[X,Y] / P2[Y,Y])
  a_res <- matrix(c(-1, beta_P), ncol=1); rownames(a_res) <- c(X,Y)
  out_res <- calc_resp_2trait(G2, P2, a_res)
  resid_ch4 <- as.numeric(out_res$R[1,1])
  resid_mbw <- as.numeric(out_res$R[2,1])
  
  # --- 3-trait G3 for ADG response (assembled) ---
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
  
  # ADG response across grid
  G_HI <- G3[, c(X,Y), drop=FALSE]
  RH <- calc_resp_H_from_B(G_HI, out$B, out$denom)
  adg <- as.numeric(RH[3,])
  
  # ADG at special points (ratio/residual)
  RH_ratio <- calc_resp_H_from_B(G_HI, out_ratio$B, out_ratio$denom)
  ratio_adg <- as.numeric(RH_ratio[3,1])
  
  RH_res <- calc_resp_H_from_B(G_HI, out_res$B, out_res$denom)
  resid_adg <- as.numeric(RH_res[3,1])
  
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

# Output 1: Monte Carlo uncertainty propagation details
write.csv(
  summary_table,
  file = "selection_index_P3/outputs/supplementary_monte_carlo_uncertainty.csv",
  row.names = FALSE
)

# Add strategy label to full grid
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

# Append ratio and residual special points as labelled rows
special_rows <- tibble(
  weight_CH4   = c(as.numeric(a_ratio[1,1]), -1),
  weight_MBW   = c(as.numeric(a_ratio[2,1]), beta_P_hat),
  response_CH4 = c(ratio_point$ch4,  resid_point$ch4),
  response_MBW = c(ratio_point$mbw,  resid_point$mbw),
  response_ADG = c(
    as.numeric(calc_resp_H_from_B(G_HI_hat, out_ratio_hat$B, out_ratio_hat$denom)[3,1]) * 1000,
    as.numeric(calc_resp_H_from_B(G_HI_hat, out_res_hat$B,   out_res_hat$denom)[3,1])   * 1000
  ),
  strategy = c("Ratio (Taylor linearisation)", "Residual (CH4 - beta*MBW)")
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
# 10) Optional: error bars (95% CI) for the two special points on the baseline cloud
# -----------------------------
ratio_ci <- list(
  ch4 = mc_summ(sim_res$ratio_ch4),
  mbw = mc_summ(sim_res$ratio_mbw)
)
resid_ci <- list(
  ch4 = mc_summ(sim_res$resid_ch4),
  mbw = mc_summ(sim_res$resid_mbw)
)

pt_ci <- bind_rows(
  tibble(point="Ratio",
         x=ratio_ci$mbw$mean, y=ratio_ci$ch4$mean,
         x_l=ratio_ci$mbw$lwr_2.5, x_u=ratio_ci$mbw$upr_97.5,
         y_l=ratio_ci$ch4$lwr_2.5, y_u=ratio_ci$ch4$upr_97.5),
  tibble(point="Residual",
         x=resid_ci$mbw$mean, y=resid_ci$ch4$mean,
         x_l=resid_ci$mbw$lwr_2.5, x_u=resid_ci$mbw$upr_97.5,
         y_l=resid_ci$ch4$lwr_2.5, y_u=resid_ci$ch4$upr_97.5)
)

# ============================================================
# 16) RESULTS-SUMMARY BLOCK (WITH MC SEs / CIs WHERE NEEDED)
#
# Assumes you have already run the Monte Carlo section and created:
#   - resp_all_adg, lr, ratio_point, resid_point, special_adg
#   - summary_table  (as shown in your output)
#
# If summary_table does not exist, this block will still print point estimates
# (but without MC SEs/CIs).
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})
# Alias for compatibility with older blocks
lr <- lr_hat

# Baseline special points table (point estimates, not MC)
special_adg <- tibble(
  point = c("Ratio (Taylor)", paste0("Residual (beta_P=", signif(beta_P_hat,4), ")")),
  ch4_response = c(ratio_point$ch4, resid_point$ch4),
  mbw_response = c(ratio_point$mbw, resid_point$mbw),
  adg_response = c(
    as.numeric(calc_resp_H_from_B(G_HI_hat, out_ratio_hat$B, out_ratio_hat$denom)[3,1]) * 1000,
    as.numeric(calc_resp_H_from_B(G_HI_hat, out_res_hat$B,   out_res_hat$denom)[3,1])   * 1000
  )
)
# -----------------------------
# Sanity checks
# -----------------------------
stopifnot(exists("resp_all_adg"), exists("lr"), exists("ratio_point"), exists("resid_point"), exists("special_adg"))
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
# 16.1 Core points table (Ratio vs Residual) + MC uncertainty
# -----------------------------
special_points <- special_adg %>%
  transmute(
    point,
    ch4 = ch4_response,
    mbw = mbw_response,
    adg = adg_response
  )

cat("\n=============================\n")
cat("CORE COMPARISON: Ratio vs Residual (point estimates)\n")
cat("=============================\n")
print(special_points)

# Add MC uncertainty (if available)
ratio_ch4_mc  <- get_mc("ratio_ch4")
ratio_mbw_mc  <- get_mc("ratio_mbw")
ratio_adg_mc  <- get_mc("ratio_adg")
resid_ch4_mc  <- get_mc("resid_ch4")
resid_mbw_mc  <- get_mc("resid_mbw")
resid_adg_mc  <- get_mc("resid_adg")

special_points_mc <- tibble(
  point = c("Ratio (Decrease CH4/MBW)", "Residual (X - beta*Y)"),
  ch4 = c(
    fmt_ci(ratio_ch4_mc$mean, ratio_ch4_mc$mc_sd, ratio_ch4_mc$lwr_2.5, ratio_ch4_mc$upr_97.5),
    fmt_ci(resid_ch4_mc$mean, resid_ch4_mc$mc_sd, resid_ch4_mc$lwr_2.5, resid_ch4_mc$upr_97.5)
  ),
  mbw = c(
    fmt_ci(ratio_mbw_mc$mean, ratio_mbw_mc$mc_sd, ratio_mbw_mc$lwr_2.5, ratio_mbw_mc$upr_97.5),
    fmt_ci(resid_mbw_mc$mean, resid_mbw_mc$mc_sd, resid_mbw_mc$lwr_2.5, resid_mbw_mc$upr_97.5)
  ),
  adg = c(
    fmt_ci(ratio_adg_mc$mean, ratio_adg_mc$mc_sd, ratio_adg_mc$lwr_2.5, ratio_adg_mc$upr_97.5),
    fmt_ci(resid_adg_mc$mean, resid_adg_mc$mc_sd, resid_adg_mc$lwr_2.5, resid_adg_mc$upr_97.5)
  ),
  n = c(ratio_ch4_mc$n, resid_ch4_mc$n)
)

cat("\n=============================\n")
cat("CORE COMPARISON: Ratio vs Residual (MC mean + uncertainty)\n")
cat("=============================\n")
if(has_mc) {
  print(special_points_mc)
} else {
  cat("summary_table not found -> MC uncertainty not printed.\n")
}

# -----------------------------
# 16.2 Feasible ADG range in the favourable quadrant (ΔMBW>0, ΔCH4<0)
# point estimates + MC uncertainty for the summary stats
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
# point estimates + MC uncertainty for the methane-optimum-with-growth metric
# -----------------------------
lr_growth_ok <- lr %>% filter(adg >= 0)

cat("\n=============================\n")
cat("GROWTH-PROTECTED REGION (ADG >= 0 within favourable quadrant)\n")
cat("=============================\n")
cat("Count:", nrow(lr_growth_ok), "out of", nrow(lr), "\n")

if(nrow(lr_growth_ok) > 0){
  best_methane_while_growth <- lr_growth_ok %>%
    arrange(ch4) %>% # most negative CH4 first (best methane reduction)
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
# 16.3 Excluded
# -----------------------------


# -----------------------------
# 16.5 Dominance checks (point-estimate geometry only)
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

ratio_row <- tibble(ch4=ratio_point$ch4, mbw=ratio_point$mbw, adg=special_adg$adg_response[1])
resid_row <- tibble(ch4=resid_point$ch4, mbw=resid_point$mbw, adg=special_adg$adg_response[2])

dom_ratio <- dominates(lr, ratio_row)
dom_resid <- dominates(lr, resid_row)

cat("\n=============================\n")
cat("DOMINANCE CHECKS (within favourable quadrant)\n")
cat("=============================\n")
cat("Points that dominate RATIO:", nrow(dom_ratio), "\n")
if(nrow(dom_ratio) > 0){
  cat("Example dominator of ratio (best ADG, then MBW, then methane):\n")
  print(dom_ratio %>% arrange(desc(adg), desc(mbw), ch4) %>% slice(1) %>%
          select(ch4, mbw, adg, w_CH4, w_MBW))
} else {
  cat("No linear-index point strictly dominates the ratio point on (CH4, MBW, ADG) simultaneously.\n")
}

cat("\nPoints that dominate RESIDUAL:", nrow(dom_resid), "\n")
if(nrow(dom_resid) > 0){
  cat("Example dominator of residual (best ADG, then MBW, then methane):\n")
  print(dom_resid %>% arrange(desc(adg), desc(mbw), ch4) %>% slice(1) %>%
          select(ch4, mbw, adg, w_CH4, w_MBW))
} else {
  cat("No linear-index point strictly dominates the residual point on (CH4, MBW, ADG) simultaneously.\n")
}

# -----------------------------
# 16.6 Index flexibility advantage (point estimate + MC range already printed)
# -----------------------------
cat("\n=============================\n")
cat("INDEX FLEXIBILITY: ADG span vs fixed points\n")
cat("=============================\n")
cat("ADG range in favourable quadrant: [",
    signif(adg_range_lr$adg_min, 4), ", ", signif(adg_range_lr$adg_max, 4), "]\n", sep="")
cat("ADG at ratio:   ", signif(special_adg$adg_response[1], 4), "\n", sep="")
cat("ADG at residual:", signif(special_adg$adg_response[2], 4), "\n", sep="")


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

# Optional: MC compact table (if you want a publication-style table of means + CIs)
if(has_mc){
  mc_compact <- tibble(
    metric = c("ratio_ch4","ratio_mbw","ratio_adg","resid_ch4","resid_mbw","resid_adg",
               "adg_min_LR","adg_p05_LR","adg_p50_LR","adg_p95_LR","adg_max_LR",
               "best_ch4_with_ADGge0","best_mbw_at_bestCH4_with_ADGge0","best_adg_at_bestCH4_with_ADGge0")
  ) %>%
    rowwise() %>%
    mutate(
      mean = get_mc(metric)$mean,
      mc_sd = get_mc(metric)$mc_sd,
      lwr_2.5 = get_mc(metric)$lwr_2.5,
      upr_97.5 = get_mc(metric)$upr_97.5,
      n = get_mc(metric)$n
    ) %>% ungroup()
  
  cat("\n=============================\n")
  cat("MONTE CARLO SUMMARY TABLE (for reporting)\n")
  cat("=============================\n")
  print(mc_compact)
} else {
  cat("\n(MC summary_table not found -> skipping MC compact table.)\n")
}








# ============================================================
# FRONTIER PLOT v3 — HULL LINE COLOURED BY ADG RESPONSE
# Key change from v2:
#   - The frontier hull line is coloured by ADG correlated response
#     rather than drawn in flat black. Each segment of the hull
#     takes the colour of the ADG value at that point.
#   - The ADG ≈ 0 dot contour layer is removed — the coloured line
#     does this job more elegantly.
#   - Background cloud points kept but slightly more transparent
#     so the coloured hull line is the dominant visual element.
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

# ── ASSUMED OBJECTS ALREADY IN ENVIRONMENT ───────────────────────────────────
# resp_all_adg, hull_df_closed, ratio_point, resid_point,
# beta_P_hat, rg_XZ_hat, rg_YZ_hat

# ── DERIVE OBJECTS ────────────────────────────────────────────────────────────
lr <- resp_all_adg %>% filter(mbw > 0, ch4 < 0)

# ── COLOUR THE HULL LINE BY ADG ───────────────────────────────────────────────
hull_lr_adg <- hull_df %>%
  filter(mbw > 0, ch4 < 0) %>%
  mutate(
    mbw_r = round(mbw, 6),
    ch4_r = round(ch4, 6)
  ) %>%
  left_join(
    resp_all_adg %>%
      mutate(mbw_r = round(mbw, 6), ch4_r = round(ch4, 6)) %>%
      select(mbw_r, ch4_r, adg),
    by = c("mbw_r", "ch4_r")
  )

# Fill any floating point misses via nearest neighbour
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

# ── SUMMARY STATISTICS ───────────────────────────────────────────────────────
adg_lims        <- range(lr$adg, na.rm = TRUE)
adg_pad         <- diff(adg_lims) * 0.05
adg_lims_padded <- c(adg_lims[1] - adg_pad, adg_lims[2] + adg_pad)
pct_positive    <- round(100 * mean(lr$adg >= 0, na.rm = TRUE), 1)
med_adg         <- median(lr$adg, na.rm = TRUE)

# ── SHARED COLOUR SCALE ───────────────────────────────────────────────────────
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

# Evenly spaced breaks from min to max, displayed to 2 dp
adg_breaks <- round(seq(adg_lims[1], adg_lims[2], length.out = 6), 2)

adg_col_scale <- scale_colour_gradientn(
  colours = palette_stops,
  values  = value_stops,
  limits  = adg_lims,
  breaks  = c(adg_lims[1], adg_lims[2]),
  labels  = sprintf("%.2f", c(adg_lims[1], adg_lims[2])),
  name    = expression("ADG correlated response (g day"^{-1}*")")
)

adg_fill_scale <- scale_fill_gradientn(
  colours = palette_stops,
  values  = value_stops,
  limits  = adg_lims,
  breaks  = c(adg_lims[1], adg_lims[2]),
  labels  = sprintf("%.2f", c(adg_lims[1], adg_lims[2])),
  name    = expression("ADG correlated response (g day"^{-1}*")")
)
# ── FIGURE ───────────────────────────────────────────────────────────────────
p_frontier <- ggplot() +
  
  # Layer 1: Background cloud
  geom_point(
    data  = lr,
    aes(x = mbw, y = ch4, colour = adg),
    alpha = 0.25,
    size  = 1.4
  ) +
  adg_col_scale +
  
  # Layer 2: Hull line coloured by ADG
  geom_path(
    data      = hull_lr_adg,
    aes(x = mbw, y = ch4, colour = adg),
    linewidth = 2.2,
    lineend   = "round",
    linejoin  = "round"
  ) +
  
  # Layer 3: Ratio special point - marker
  geom_point(
    data   = ratio_point,
    aes(x = mbw, y = ch4),
    size   = 5,
    shape  = 24,
    fill   = "white",
    colour = "#1A1A2E",
    stroke = 1.5
  ) +
  annotate(
    "richtext",
    x        = ratio_point$mbw + 0.02,
    y        = ratio_point$ch4 - 0.05,
    label    = paste0(
      "<b>Ratio</b><br>",
      "\u0394CH4: ", round(ratio_point$ch4, 2), " (\u00b1", round(get_mc("ratio_ch4")$mc_sd, 2), ")<br>",
      "\u0394MBW: ", round(ratio_point$mbw, 2),  " (\u00b1", round(get_mc("ratio_mbw")$mc_sd, 2), ")<br>",
      "\u0394ADG: ", round(special_adg$adg_response[1], 2), " (\u00b1", round(get_mc("ratio_adg")$mc_sd * 1000, 2), ")"
    ),
    size      = 3.3,
    colour    = "#1A1A2E",
    hjust     = 0,
    fill      = NA,
    label.color = NA
  ) +
  
  # Layer 4: Residual special point - marker
  geom_point(
    data   = resid_point %>% filter(mbw > 0, ch4 < 0),
    aes(x = mbw, y = ch4),
    size   = 5,
    shape  = 22,
    fill   = "white",
    colour = "#1A1A2E",
    stroke = 1.5
  ) +
  annotate(
    "richtext",
    x        = resid_point$mbw + 0.02,
    y        = resid_point$ch4 - 0.05,
    label    = paste0(
      "<b>Residual</b><br>",
      "\u0394CH4: ", round(resid_point$ch4, 2), " (\u00b1", round(get_mc("resid_ch4")$mc_sd, 2), ")<br>",
      "\u0394MBW: ", round(resid_point$mbw, 2),  " (\u00b1", round(get_mc("resid_mbw")$mc_sd, 2), ")<br>",
      "\u0394ADG: ", round(special_adg$adg_response[2], 2), " (\u00b1", round(get_mc("resid_adg")$mc_sd * 1000, 2), ")"
    ),
    size      = 3.3,
    colour    = "#1A1A2E",
    hjust     = 0,
    fill      = NA,
    label.color = NA
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
  "selection_index_P3/outputs/frontier_adg_v5.png",
  width  = 8,
  height = 8,
  dpi    = 300,
  bg     = "white"
)




