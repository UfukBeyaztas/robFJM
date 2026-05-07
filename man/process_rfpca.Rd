\name{process_rfpca}
\alias{process_rfpca}
\title{
Process Functional Predictors by Robust FPCA
}
\description{
Computes robust functional principal component scores and basis evaluations
for functional predictors observed on a common grid.
}
\usage{
process_rfpca(subjs, funcs_mat, gp, nbasis = 15, npc_max = NULL)
}
\arguments{
  \item{subjs}{
A vector of subject identifiers. Its length must match the number of rows of
\code{funcs_mat}. If repeated rows are present for the same subject, robust
FPCA is computed using the first occurrence of each subject and the resulting
scores are mapped back to all rows according to \code{subjs}.
}
  \item{funcs_mat}{
A numeric matrix containing the observed functional predictor values. Rows
correspond to observations or subjects and columns correspond to grid points.
The number of columns must match the length of \code{gp}.
}
  \item{gp}{
A numeric vector giving the grid points on which the functional predictor is
observed.
}
  \item{nbasis}{
An integer specifying the number of B-spline basis functions used to smooth and
represent the functional observations before robust FPCA. The default is
\code{15}.
}
  \item{npc_max}{
An optional integer specifying the maximum number of robust functional
principal components to retain. If \code{NULL}, the number of retained
components is selected automatically by the robust FPCA step.
}
}
\details{
This function is the core robust functional dimension-reduction routine used
by the proposed robust functional joint modeling procedure. It first identifies
the subject-level functional observations by retaining the first occurrence of
each subject identifier in \code{subjs}. Robust FPCA is then applied to these
baseline subject-level functional observations through the internal routine
\code{getRPCA}.

The robust FPCA step represents each function in a B-spline basis and computes
robust principal directions using a projection-pursuit covariance estimator.
Let \eqn{\mathcal{X}_i(s)} denote a functional predictor observed on the grid
\code{gp}. The method approximates
\deqn{
  \mathcal{X}_i(s) \approx \sum_{\ell=1}^{K_b} c_{i\ell}\psi_\ell(s),
}
where \eqn{\{\psi_\ell\}_{\ell=1}^{K_b}} is a B-spline basis. A robust center
is used for the spline coefficient vectors, and robust principal directions are
then extracted from a projection-pursuit covariance estimate. This approach is
designed to reduce the influence of atypical curves and local functional
artifacts on the estimated eigenfunctions.

If \code{npc_max = NULL}, the function uses the number of robust components
selected by \code{getRPCA}. If \code{npc_max} is supplied, the number of
retained components is
\deqn{
  K = \min(K_{\mathrm{rob}}, \code{npc\_max}),
}
where \eqn{K_{\mathrm{rob}}} is the number of components selected by the robust
FPCA step. At least one component is retained.

The output score matrix \code{FPCscore} has the same number of rows as
\code{funcs_mat}. Therefore, when the input contains repeated longitudinal
records for the same subject, each repeated record receives the same
subject-level robust FPCA score vector. This convention is useful in the
functional joint model because the functional predictor is treated as a
baseline subject-level predictor, while the longitudinal scalar response may
be measured repeatedly over time.

The returned basis evaluations in \code{basis} are evaluated on \code{gp} and
are used by \code{\link{fit_FAPHL}} to reconstruct the estimated functional
coefficient functions after score-level model fitting.
}
\value{
A list with the following components:
  \item{FPCscore}{
  A numeric matrix of robust functional principal component scores. The number
  of rows equals the number of rows of \code{funcs_mat}, and the number of
  columns equals the number of retained robust components. Columns are named
  \code{PC1}, \code{PC2}, and so on.
  }
  \item{kz}{
  The number of retained robust functional principal components.
  }
  \item{basis}{
  A numeric matrix containing the retained robust eigenfunctions evaluated on
  the grid \code{gp}. Rows correspond to grid points and columns correspond to
  retained components.
  }
  \item{rp}{
  The full robust FPCA object returned by \code{getRPCA}. This object contains
  additional information such as the functional principal component object,
  robust scores, B-spline basis, and grid information.
  }
}
\note{
This function is usually not called directly by end users. The recommended
workflow is to call \code{\link{add_scores}}, which applies
\code{process_rfpca} separately to the longitudinal-side and survival-side
functional predictors and attaches the resulting robust scores to the data.

The function assumes that the functional predictor is observed on a common
grid for all subjects. Missing values in \code{funcs_mat} should be handled
before calling this function.
}
\seealso{
\code{\link{add_scores}}, \code{\link{fit_FAPHL}},
\code{\link{simulate_data}}
}
\examples{
\donttest{

I <- 100; J <- 4; D <- 10; by <- 0.1

# Contamination settings
outlier_frac_Y   <- 0.05
outlier_frac_fun <- 0.05

dat <- simulate_data(
  I = I, J = J, D = D, by = by,
  outlier_frac_Y   = outlier_frac_Y,
  outlier_frac_fun = outlier_frac_fun
)

sgrid <- seq(0, D, by)
m <- length(sgrid)
grid_cols_X <- paste0("func.X.", 1:m)

Xmat <- as.matrix(dat$data.long[, grid_cols_X, drop = FALSE])

rfpca_x <- process_rfpca(
  subjs = dat$data.long$ID,
  funcs_mat = Xmat,
  gp = dat$grid,
  nbasis = 10
)

names(rfpca_x)
dim(rfpca_x$FPCscore)
rfpca_x$kz

matplot(dat$grid, rfpca_x$basis, type = "l",
        xlab = "Functional domain",
        ylab = "Robust eigenfunctions")

}
}
