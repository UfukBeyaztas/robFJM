\name{fit_FAPHL}
\alias{fit_FAPHL}
\title{
Fit a Robust Functional Joint Model
}
\description{
Fits the proposed robust functional joint model using robust FPCA scores and
Student-\eqn{t} adjusted profile \eqn{h}-likelihood estimation.
}
\usage{
fit_FAPHL(long, surv, FPCA.X, FPCA.W, by,
          nu = 4, verbose = FALSE)
}
\arguments{
  \item{long}{
A longitudinal data frame containing repeated scalar responses and robust FPCA
scores for the longitudinal-side functional predictor. The data frame must
include the columns \code{ID}, \code{Y}, \code{obstime}, and \code{FPCscore.X}.
The column \code{ID} is the subject identifier, \code{Y} is the scalar
longitudinal response, \code{obstime} is the longitudinal observation time, and
\code{FPCscore.X} is the robust FPCA score representation of the functional
predictor entering the longitudinal submodel.
}
  \item{surv}{
A subject-level survival data frame. The data frame must include the columns
\code{ID}, \code{time}, \code{event}, \code{W1}, and \code{FPCscore.W}.
The column \code{time} is the observed event or censoring time, \code{event}
is the event indicator, \code{W1} is a scalar survival-side explanatory
variable, and \code{FPCscore.W} is the robust FPCA score representation of the
functional predictor entering the survival submodel.
}
  \item{FPCA.X}{
The robust FPCA object for the longitudinal-side functional predictor. This is
typically obtained from \code{\link{add_scores}} and should contain the robust
basis evaluations in \code{FPCA.X$basis}.
}
  \item{FPCA.W}{
The robust FPCA object for the survival-side functional predictor. This is
typically obtained from \code{\link{add_scores}} and should contain the robust
basis evaluations in \code{FPCA.W$basis}.
}
  \item{by}{
A numeric value giving the spacing of the functional grid. This argument is
kept for consistency with the simulation and score-construction interface. The
coefficient functions are reconstructed on the grid represented in
\code{FPCA.X$basis} and \code{FPCA.W$basis}.
}
  \item{nu}{
Degrees of freedom parameter for the Student-\eqn{t} longitudinal error
distribution. Smaller values yield heavier tails and stronger downweighting of
large longitudinal residuals. The default is \code{4}.
}
  \item{verbose}{
Logical argument reserved for optional progress output. The default is
\code{FALSE}.
}
}
\details{
The function implements the main estimation step of the proposed robust
functional joint model. It assumes that robust functional dimension reduction
has already been performed, usually by \code{\link{add_scores}}. The
longitudinal and survival functional predictors are therefore represented by
robust FPCA scores, which enter the two model components as finite-dimensional
explanatory variables.

The fitted longitudinal submodel has the working form
\deqn{
  Y_{ij}
  =
  \beta_0 + \beta_1 t_{ij}
  + \boldsymbol{\xi}_{i}^{(\mathcal{X})\top}\boldsymbol{\beta}_x
  + u_i + \varepsilon_{ij},
}
where \eqn{Y_{ij}} is the scalar longitudinal response, \eqn{t_{ij}} is the
observation time, \eqn{\boldsymbol{\xi}_{i}^{(\mathcal{X})}} is the robust FPCA
score vector for the longitudinal-side functional predictor, and \eqn{u_i} is
a shared subject-specific random effect.

The survival submodel uses a Weibull proportional hazards structure with a
current-value association term,
\deqn{
  h_i(t)
  =
  \tau t^{\tau - 1}
  \exp\{
    \gamma_1 W_{1i}
    + \boldsymbol{\xi}_{i}^{(\mathcal{W})\top}\boldsymbol{\gamma}_w
    + \alpha m_i(t)
  \},
}
where \eqn{\boldsymbol{\xi}_{i}^{(\mathcal{W})}} is the robust FPCA score
vector for the survival-side functional predictor, \eqn{m_i(t)} is the latent
longitudinal trajectory, \eqn{\alpha} is the association parameter, and
\eqn{\tau} is the Weibull shape parameter.

Robustness to aberrant longitudinal observations is introduced through a
Student-\eqn{t} error distribution. Using the Gaussian scale-mixture
representation,
\deqn{
  \varepsilon_{ij} \mid \lambda_{ij}
  \sim
  N\left(0, \sigma_\varepsilon^2/\lambda_{ij}\right),
  \qquad
  \lambda_{ij}
  \sim
  \mathrm{Gamma}(\nu/2,\nu/2),
}
the method updates observation-specific weights according to
\deqn{
  \widehat{\lambda}_{ij}
  =
  \frac{\nu + 1}{\nu + r_{ij}^2/\sigma_\varepsilon^2},
}
where \eqn{r_{ij}} is the current residual. Large residuals therefore receive
smaller weights.

The function then calls the internal routine \code{fit_APHL}, which maximizes
an adjusted profile \eqn{h}-likelihood. The random effects are profiled by
damped Newton updates, survival integrals are evaluated numerically, and the
Laplace-type curvature correction is used to adjust the profiled
\eqn{h}-likelihood. After estimating the score-level coefficients, the
function maps them back to the functional domain to obtain
\eqn{\widehat{\beta}(s)} and \eqn{\widehat{\gamma}(s)}.
}
\value{
A list with the following components:
  \item{fit}{
  A list returned by the internal adjusted profile \eqn{h}-likelihood fitting
  routine. It contains the estimated standardized-model parameters after
  back-transformation, including fixed-effect vectors, Weibull shape parameter,
  association parameter, variance components, and estimated random effects.
  }
  \item{BetaHat}{
  A numeric vector containing the reconstructed estimate of the functional
  coefficient \eqn{\beta(s)} for the longitudinal-side functional predictor,
  evaluated on the grid represented by \code{FPCA.X$basis}.
  }
  \item{GammaHat}{
  A numeric vector containing the reconstructed estimate of the functional
  coefficient \eqn{\gamma(s)} for the survival-side functional predictor,
  evaluated on the grid represented by \code{FPCA.W$basis}.
  }
  \item{scalars}{
  A named numeric vector containing the main scalar parameter estimates:
  \code{beta1}, \code{var_b}, \code{var_eps}, \code{gamma1}, \code{alpha}, and
  \code{Tau}. Here \code{beta1} is the longitudinal time effect, \code{var_b}
  is the random-intercept variance estimate, \code{var_eps} is the
  longitudinal error variance estimate, \code{gamma1} is the scalar
  survival-side coefficient for \code{W1}, \code{alpha} is the current-value
  association parameter, and \code{Tau} is the Weibull shape parameter.
  }
}
\note{
This function expects the robust FPCA score columns to be present in the input
data. The recommended workflow is to first call \code{\link{simulate_data}} or
prepare user data in the required format, then call \code{\link{add_scores}},
and finally fit the proposed model with \code{fit_FAPHL}.

The current implementation uses a random-intercept shared-effects structure and
a Weibull baseline hazard. The scalar survival explanatory variable is expected
to be named \code{W1}.
}
\seealso{
\code{\link{add_scores}}, \code{\link{process_rfpca}},
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

fit <- fit_FAPHL(
  long = sc$robust$long,
  surv = sc$robust$surv,
  FPCA.X = sc$robust$FPCA.X,
  FPCA.W = sc$robust$FPCA.W,
  by = by, verbose = FALSE
)

fit$scalars
length(fit$BetaHat)
length(fit$GammaHat)

plot(dat$grid, fit$BetaHat, type = "l",
     xlab = "Functional domain",
     ylab = "Estimated beta(s)")

plot(dat$grid, fit$GammaHat, type = "l",
     xlab = "Functional domain",
     ylab = "Estimated gamma(s)")
}
}
