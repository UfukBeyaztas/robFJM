# robFJM <img src="https://img.shields.io/badge/R-%3E=3.5.0-1f425f.svg" alt="R (>= 3.5.0)" align="right" height="20"/>

**robFJM** provides tools for fitting **robust functional joint models** for data involving **longitudinal scalar responses**, **time-to-event outcomes**, and **functional predictors**. The package is designed for studies in which repeated scalar outcomes and event times are jointly modeled together with densely observed or high-dimensional functional predictors.

The package implements a two-stage robust estimation framework. First, functional predictors entering the longitudinal and survival submodels are represented by **robust functional principal component analysis** (RFPCA), reducing the influence of atypical curves and functional artifacts. Second, conditional on these robust score representations, the joint model is fitted using **Student-\(t\) adjusted profile \(h\)-likelihood** estimation with observation-specific downweighting, shared random effects, a Weibull baseline hazard, numerical quadrature for survival integrals, and reconstruction of the longitudinal and survival functional coefficient functions.

The package also includes simulation tools for generating clean and contaminated functional joint model data, allowing users to evaluate the robustness of the proposed method under scalar-response outliers and functional-predictor contamination.

---

## 🚀 Key Features

- **Robust functional joint modeling:** fits joint models for scalar longitudinal outcomes, survival outcomes, and functional predictors.

- **Two functional predictors:** allows separate functional predictors to enter the longitudinal and survival submodels.

- **Robust functional dimension reduction:** uses robust FPCA to obtain stable score representations for functional predictors.

- **Resistance to functional contamination:** reduces the influence of atypical curves, magnitude-shifted functions, and local functional artifacts.

- **Student-\(t\) longitudinal errors:** downweights aberrant scalar longitudinal responses through observation-specific robust weights.

- **Adjusted profile \(h\)-likelihood estimation:** profiles shared random effects and applies curvature correction for likelihood-based frequentist estimation.

- **Weibull survival component:** fits a Weibull proportional hazards structure with a current-value association term linking the longitudinal trajectory to event risk.

- **Shared random effects:** accounts for subject-specific dependence in repeated longitudinal measurements.

- **Functional coefficient reconstruction:** reconstructs the longitudinal coefficient function \(\widehat{\beta}(s)\) and survival coefficient function \(\widehat{\gamma}(s)\) from estimated robust score-level coefficients.

- **Built-in simulation tools:** generates clean and contaminated data from a functional joint model with longitudinal, survival, and functional components.

- **User manual included:** the package is accompanied by the manual file `robFJM_1.0.0.pdf`.

---

## 📦 Installation

You can install the development version of **robFJM** from GitHub:

```r
install.packages("remotes")
remotes::install_github("UfukBeyaztas/robFJM")

Then load the package:

library(robFJM)
📘 Main Functions

The package contains four main user-facing functions:

Function	Description
simulate_data()	Simulates longitudinal scalar responses, survival outcomes, and functional predictors under clean or contaminated functional joint model settings.
add_scores()	Computes robust FPCA scores for the longitudinal-side and survival-side functional predictors and attaches them to the data.
process_rfpca()	Performs robust FPCA for functional predictors observed on a common grid.
fit_FAPHL()	Fits the proposed robust functional joint model using robust FPCA scores and Student-t adjusted profile h-likelihood estimation.
