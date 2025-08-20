\name{Filtros}
\alias{get_annual_cycle}
\alias{filter_annual_cycle}
\alias{fit_all_coef}
\alias{filter_sin}
\alias{fit_all_coef_with_sig}
\alias{filter_sin_sigmoide1}
\alias{filter_sin_sigmoide2}
\title{Muchos Filtros}
\usage{
get_annual_cycle(x,date)
filter_annual_cycle(t,annual_cycle)
fit_all_coef(x,date, k_sim = 4)
filter_sin(t, coef)
fit_all_coef_with_sig(x,date, k_sim_sig = 2)
filter_sin_sigmoide1(t, coef)
filter_sin_sigmoide2(t, coef)
}
\arguments{
  \item{x}{Vector with the Julian Day associated to the PP_vec}
  \item{date}{Numeric Vector for which percentile will be compute}
  \item{t}{numeric, percentile to be computed}
  \item{annual_cycle}{integer, window size}
  \item{coef}{integer, window size}
  \item{k_sim}{integer, window size}
  \item{k_sim_sig}{integer, window size}

}
\description{
asdasd
}

