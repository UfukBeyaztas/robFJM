\name{simulate_data}
\alias{simulate_data}
\title{
Simulate Functional Joint Model Data
}
\description{
Generates synthetic longitudinal, survival, and functional predictor data for
illustrating and evaluating the proposed robust functional joint model.
}
\usage{
simulate_data(I = 200, J = 4,
              D = 10, by = 0.1,
              beta0 = 1.35, beta1 = 0.78,
              var_b = 1, var_eps = 1,
              gamma0 = 0, gamma1 = -1.75, alpha = 0.29,
              Tau = 0.8,
              lambdaC = -6.0,
              Kx_latent = 5, Kw_latent = 5,
              meas_sd_X = 0.15, meas_sd_W = 0.15,
              outlier_frac_Y = 0.00, outlier_mag_Y = c(20, 30),
              outlier_frac_fun = 0.00, outlier_mag_fun = c(3, 6))
}
\arguments{
  \item{I}{
An integer specifying the number of subjects to generate. The default is
\code{200}.
}
  \item{J}{
An integer specifying the planned number of repeated longitudinal measurements
per subject before trimming measurements after the observed follow-up time. The
default is \code{4}.
}
  \item{D}{
A numeric value giving the upper limit of the functional domain and the planned
longitudinal follow-up window. The functional grid is defined on
\eqn{[0,D]}. The default is \code{10}.
}
  \item{by}{
A numeric value giving the spacing between adjacent functional grid points. The
grid is generated as \code{seq(0, D, by)}. The default is \code{0.1}.
}
  \item{beta0}{
The intercept parameter in the longitudinal submodel. The default is
\code{1.35}.
}
  \item{beta1}{
The longitudinal time-effect parameter. The default is \code{0.78}.
}
  \item{var_b}{
The variance of the subject-specific random intercept in the longitudinal
submodel. The default is \code{1}.
}
  \item{var_eps}{
The variance of the scalar longitudinal measurement error. The default is
\code{1}.
}
  \item{gamma0}{
The intercept term used in the survival linear predictor. The default is
\code{0}.
}
  \item{gamma1}{
The scalar survival-side coefficient for the binary explanatory variable
\code{W1}. The default is \code{-1.75}.
}
  \item{alpha}{
The current-value association parameter linking the latent longitudinal
trajectory to the hazard function. The default is \code{0.29}.
}
  \item{Tau}{
The Weibull shape parameter used in both event-time and censoring-time
generation. The default is \code{0.8}.
}
  \item{lambdaC}{
A log-scale censoring parameter controlling the censoring-time distribution.
Larger values generally produce shorter censoring times and hence more
censoring. The default is \code{-6.0}.
}
  \item{Kx_latent}{
The number of latent basis functions used to generate the functional predictor
entering the longitudinal submodel. The default is \code{5}.
}
  \item{Kw_latent}{
The number of latent basis functions used to generate the functional predictor
entering the survival submodel. The default is \code{5}.
}
  \item{meas_sd_X}{
The standard deviation of the measurement noise added to the longitudinal-side
functional predictor. The default is \code{0.15}.
}
  \item{meas_sd_W}{
The standard deviation of the measurement noise added to the survival-side
functional predictor. The default is \code{0.15}.
}
  \item{outlier_frac_Y}{
The proportion of scalar longitudinal observations contaminated by additive
outliers. The default is \code{0.00}, corresponding to no scalar-response
contamination.
}
  \item{outlier_mag_Y}{
A numeric vector of length two giving the lower and upper bounds for the
absolute magnitude of scalar longitudinal outliers. The sign is randomly
assigned. The default is \code{c(20, 30)}.
}
  \item{outlier_frac_fun}{
The proportion of subjects whose functional predictors are contaminated by
subject-level magnitude shifts. The default is \code{0.00}, corresponding to no
functional contamination.
}
  \item{outlier_mag_fun}{
A numeric vector of length two giving the lower and upper bounds for the
absolute magnitude of functional-predictor shifts. The sign is randomly
assigned. The default is \code{c(3, 6)}.
}
}
\details{
This function generates data from a functional joint model with a scalar
longitudinal response, a time-to-event outcome, and two baseline functional
predictors. One functional predictor enters the longitudinal submodel and the
other enters the survival submodel.

The longitudinal submodel is generated according to
\deqn{
  Y_{ij}
  =
  \beta_0 + \beta_1 t_{ij}
  + \int_0^D X_i(s)\beta(s)\,ds
  + b_i + \varepsilon_{ij},
}
where \eqn{b_i \sim N(0,\sigma_b^2)} is a subject-specific random intercept
and \eqn{\varepsilon_{ij} \sim N(0,\sigma_\varepsilon^2)} is the scalar
longitudinal error. In the implementation, \eqn{\sigma_b^2} is controlled by
\code{var_b} and \eqn{\sigma_\varepsilon^2} is controlled by \code{var_eps}.

The survival submodel uses a Weibull current-value joint modeling structure.
The subject-specific hazard is generated from
\deqn{
  h_i(t)
  =
  \tau t^{\tau-1}
  \exp\{
    \gamma_0 + \gamma_1 W_{1i}
    + \int_0^D W_i(s)\gamma(s)\,ds
    + \alpha m_i(t)
  \},
}
where \eqn{m_i(t)} is the latent longitudinal trajectory. Event times are
generated by numerical inversion of the cumulative hazard. Censoring times are
generated independently from a Weibull-type distribution controlled by
\code{lambdaC} and \code{Tau}.

The true coefficient functions are
\deqn{
  \beta(s) = 2\sin(\pi s/5),
  \qquad
  \gamma(s) = 1.2\sin(\pi s/4).
}
The latent functional predictors are generated using finite sine and cosine
basis expansions with decreasing score variances. Measurement noise is then
added to obtain the observed functional predictors.

Two contamination mechanisms can be introduced. First, scalar longitudinal
responses can be contaminated by additive outliers using \code{outlier_frac_Y}
and \code{outlier_mag_Y}. Second, functional predictors can be contaminated by
subject-level vertical shifts using \code{outlier_frac_fun} and
\code{outlier_mag_fun}. These mechanisms are useful for evaluating the
robustness of the proposed robust FPCA plus Student-\eqn{t} adjusted profile
\eqn{h}-likelihood estimator.

After generating the full longitudinal data, observations occurring after the
observed event or censoring time are removed. Thus, the returned longitudinal
data mimic the usual joint-modeling setting in which repeated measurements are
only available up to follow-up.
}
\value{
A list with the following components:
\describe{
  \item{data.long}{
  A longitudinal data frame containing repeated measurements. It includes
  subject identifiers \code{ID}, visit numbers \code{visit}, observation times
  \code{obstime}, observed follow-up times \code{time}, event indicators
  \code{event}, scalar longitudinal responses \code{Y}, functional predictor
  columns named \code{func.X.1}, \code{func.X.2}, and so on, and the generated
  random intercept \code{b}.
  }
  \item{data.surv}{
  A subject-level survival data frame containing \code{ID}, observed follow-up
  time \code{time}, event indicator \code{event}, scalar explanatory variable
  \code{W1}, and functional predictor columns named \code{func.W.1},
  \code{func.W.2}, and so on.
  }
  \item{grid}{
  The numeric grid on which the functional predictors and true coefficient
  functions are evaluated.
  }
  \item{by}{
  The functional grid spacing used to generate \code{grid}.
  }
  \item{true}{
  A list containing the true parameter values and true coefficient functions
  used in the data-generating process. It includes \code{beta0}, \code{beta1},
  \code{var_b}, \code{var_eps}, \code{gamma0}, \code{gamma1}, \code{alpha},
  \code{Tau}, \code{trueBeta}, and \code{trueGamma}.
  }
}
}
\seealso{
\code{\link{add_scores}}, \code{\link{process_rfpca}},
\code{\link{fit_FAPHL}}
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

names(dat)
head(dat$data.long)
head(dat$data.surv)
}
}
