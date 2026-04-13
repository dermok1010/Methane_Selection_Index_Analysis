# Selection Response Comparison of Methane Traits with Monte Carlo Uncertainty

## Overview

This repository contains an R-based quantitative genetics framework for comparing **alternative methane breeding objectives** using selection index theory and Monte Carlo uncertainty propagation.

The script evaluates how different methane trait definitions influence predicted genetic responses in:

- **Methane production (CH4/day)**
- **Metabolic body weight (MBW)**
- **Average daily gain (ADG)**

It contrasts traditional derived traits such as:

- **Methane ratio traits** (e.g. CH4 / MBW)
- **Residual methane traits** (CH4 adjusted for MBW)
- **Flexible linear breeding indexes**

The aim is to show how different trait constructions implicitly constrain multi-trait selection responses, while linear indexes allow broader optimisation across production and emissions traits.

---

## Scientific Motivation

Methane efficiency traits are often expressed as ratios or residuals. However, these derived traits impose fixed weighting structures between component traits.

For example:

- Selecting on **CH4 / MBW** rewards larger body size while penalising methane.
- Selecting on **Residual CH4** targets methane independent of MBW.
- A **general linear index** allows any weighting between methane reduction and growth traits.

This script quantifies the achievable response space and demonstrates whether ratio or residual traits sit on, inside, or outside the wider multi-trait frontier.

---

## Methods Summary

## 1. Variance Component Inputs

The script uses variance components and standard errors extracted from **ASReml `.pvc` outputs**, including:

- Additive genetic variances
- Residual variances
- Permanent environmental variances
- Genetic correlations
- Residual correlations

Traits analysed:

- CH4 production
- Metabolic body weight
- Average daily gain

---

## 2. Selection Index Theory

Optimal index weights are computed using:

\[
b = P^{-1}Ga
\]

Where:

- **G** = additive genetic covariance matrix  
- **P** = phenotypic covariance matrix  
- **a** = desired breeding objective vector  
- **b** = optimal measurable index weights

Predicted correlated responses are then calculated for all traits.

---

## 3. Response Frontier

A grid of methane vs body-weight objectives is explored to generate the full achievable response cloud.

The script then identifies the **convex hull frontier**, representing Pareto-optimal trade-offs where one trait cannot improve without worsening another.

---

## 4. Special Trait Comparisons

Two common methane objectives are mapped onto the frontier:

### Ratio Objective

Approximate selection on:

\[
CH4 / MBW
\]

using first-order Taylor linearisation.

### Residual Objective

Selection on:

\[
CH4 - \beta MBW
\]

where β is the phenotypic regression coefficient.

---

## 5. Monte Carlo Uncertainty

Uncertainty in ASReml estimates is propagated using 1,000 Monte Carlo simulations.

Sampling uses:

- Log-normal draws for variances
- Fisher-z transformed draws for correlations

Outputs include:

- Means
- Monte Carlo SD
- 95% confidence intervals

---

## Key Outputs

## Figures

### Frontier Plot

Shows achievable combinations of:

- ΔCH4
- ΔMBW

with hull coloured by correlated response in ADG.

Special points:

- Triangle = Ratio trait
- Square = Residual trait

Saved as:

```r
frontier_adg_v5.png
