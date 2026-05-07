\name{add_scores}
\alias{add_scores}
\title{
Add Robust Functional Principal Component Scores
}
\description{
Adds robust functional principal component score representations to the
longitudinal and survival data sets used in the proposed robust functional
joint modeling framework.
}
\usage{
add_scores(long, surv, grid_cols_X, grid_cols_W,
           gp, rfpca_nbasis = 15, rfpca_npc_max = NULL)
}
\arguments{
  \item{long}{
A data frame containing the longitudinal observations. The data frame must
include a subject identifier column named \code{ID}, together with the
functional predictor columns specified by \code{grid_cols_X}. In the proposed
functional joint model, these functional measurements enter the longitudinal
submodel.
}
  \item{surv}{
A data frame containing the subject-level survival information. The data frame
must include a subject identifier column named \code{ID}, together with the
functional predictor columns specified by \code{grid_cols_W}. In the proposed
functional joint model, these functional measurements enter the survival
submodel.
}
  \item{grid_cols_X}{
A character vector or integer vector identifying the columns of \code{long}
that contain the longitudinal-side functional predictor values observed on the
common grid \code{gp}.
}
  \item{grid_cols_W}{
A character vector or integer vector identifying the columns of \code{surv}
that contain the survival-side functional predictor values observed on the
common grid \code{gp}.
}
  \item{gp}{
A numeric vector giving the grid points at which the functional predictors are
observed. The length of \code{gp} must match the number of functional grid
columns selected by \code{grid_cols_X} and \code{grid_cols_W}.
}
  \item{rfpca_nbasis}{
An integer specifying the number of B-spline basis functions used in the robust
functional principal component analysis. The default is \code{15}.
}
  \item{rfpca_npc_max}{
An optional integer giving the maximum number of robust functional principal
components to retain. If \code{NULL}, the number of retained components is
chosen automatically by the robust FPCA step.
}
}
\details{
The function is a preprocessing wrapper for the proposed robust functional
joint model. It extracts the functional predictor matrix from the longitudinal
data and the functional predictor matrix from the survival data, and then
applies robust functional principal component analysis separately to the two
functional predictors.

Let \eqn{\mathcal{X}_i(s)} denote the functional predictor entering the
longitudinal submodel and \eqn{\mathcal{W}_i(s)} denote the functional
predictor entering the survival submodel. The proposed method represents the
functional linear terms
\deqn{
  \int \mathcal{X}_i(s)\beta(s)\,ds
  \quad \mbox{and} \quad
  \int \mathcal{W}_i(s)\gamma(s)\,ds
}
through low-dimensional robust score vectors. Specifically, robust FPCA is
used to obtain estimated eigenfunctions
\eqn{\widehat{\phi}_1^{(\mathcal{X},R)}, \ldots,
\widehat{\phi}_{K_{\mathcal{X}}}^{(\mathcal{X},R)}} and
\eqn{\widehat{\phi}_1^{(\mathcal{W},R)}, \ldots,
\widehat{\phi}_{K_{\mathcal{W}}}^{(\mathcal{W},R)}}. The corresponding robust
scores are then attached to the data frames as list-columns named
\code{FPCscore.X} and \code{FPCscore.W}.

In contrast to classical FPCA, the robust FPCA step is designed to reduce the
effect of atypical curves and local functional artifacts on the extracted
score representation. These robust scores are then used as finite-dimensional
predictors in \code{\link{fit_FAPHL}}, where the joint model is estimated by
the proposed Student-\eqn{t} adjusted profile \eqn{h}-likelihood method.

The current implementation returns a nested list with one top-level element,
\code{robust}. Therefore, downstream calls should use objects such as
\code{out$robust$long}, \code{out$robust$surv}, \code{out$robust$FPCA.X}, and
\code{out$robust$FPCA.W}.
}
\value{
A list with the following component:
\describe{
  \item{robust}{
  A list containing the robust FPCA-transformed data and fitted robust FPCA
  objects:
  \describe{
    \item{long}{
    A copy of \code{long} with an added list-column \code{FPCscore.X}
    containing the robust FPCA scores for the longitudinal-side functional
    predictor.
    }
    \item{surv}{
    A copy of \code{surv} with an added list-column \code{FPCscore.W}
    containing the robust FPCA scores for the survival-side functional
    predictor.
    }
    \item{FPCA.X}{
    The robust FPCA object for the longitudinal-side functional predictor.
    This object contains the robust scores, retained basis evaluations, number
    of retained components, and additional robust FPCA output.
    }
    \item{FPCA.W}{
    The robust FPCA object for the survival-side functional predictor. This
    object contains the robust scores, retained basis evaluations, number of
    retained components, and additional robust FPCA output.
    }
  }
  }
}
}
\note{
This function is intended to be used before \code{\link{fit_FAPHL}}. The input
data should be organized so that each subject has a unique \code{ID}; the
longitudinal data may contain repeated rows per subject, whereas the survival
data should contain one row per subject.
}
\seealso{
\code{\link{process_rfpca}}, \code{\link{fit_FAPHL}},
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
grid_cols_W <- paste0("func.W.", 1:m)

sc <- add_scores(
  long = dat$data.long,
  surv = dat$data.surv,
  grid_cols_X = grid_cols_X,
  grid_cols_W = grid_cols_W,
  gp = dat$grid,
  rfpca_nbasis = 10
)

names(sc)
names(sc$robust)
dim(sc$robust$FPCA.X$FPCscore)
dim(sc$robust$FPCA.W$FPCscore)
}
}
