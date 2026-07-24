\name{simulate_sft}
\alias{simulate_sft}
\alias{simulate_sft_group}
\alias{find_equal_vc_yes}
\title{SFT Forward Simulation}
\description{
  Tools for simulating full Systems Factorial Technology experimental designs 
  using cognitive process models. \code{simulate_sft} simulates a single participant, 
  and \code{simulate_sft_group} generates hierarchical data for multiple participants. 
  \code{find_equal_vc_yes} is a helper for LBA capacity matching.
}
\usage{
simulate_sft(model = c("lba", "ou"), n, p_vec, design = NULL, salience = FALSE, logical_rules = NULL, sics = FALSE, cdfs = FALSE, a_t = FALSE, subject = 1, backend = c("auto", "Rcpp", "R"), dt = .001, max_time = 2, bridge = TRUE, bridge_exp_cutoff = -20, adaptive = TRUE, adaptive_mode = c("ptol", "distance"), adapt_factor = 32, adapt_eps = 1.5, p_tol = 1e-8, curv_tol = .05, kappa = NULL, tau = NULL, rho = NULL, architecture = c("parallel", "serial"), serial_mode = c("self-terminating", "exhaustive"), serial_order = c("random", "A", "B", "fixed"), p_A_first = .5, serial_first = NULL, serial_stopping = NULL, serial_p_A_first = NULL, serial_type = NULL)
simulate_sft_group(model = c("lba", "ou"), n, params, nSubjects = NULL, subject_ids = NULL, combine = TRUE, ...)
find_equal_vc_yes(vno, t, A, b, t0, sv, interval = c(1, 4.5))
}
\arguments{
  \item{model}{The process model to simulate: \code{"lba"} (Linear Ballistic Accumulator) or \code{"ou"} (Ornstein-Uhlenbeck).}
  \item{n}{Number of simulated trials per design cell.}
  \item{p_vec}{Named parameter vector defining the simulated process.}
  \item{design}{Design specifications for factors and levels.}
  \item{salience}{Logical; whether to include salience manipulations for SIC tests.}
  \item{logical_rules}{Rules (\code{"OR"}, \code{"AND"}) to simulate over.}
  \item{sics, cdfs, a_t}{Output toggles for automatically computing theoretical curves.}
  \item{subject}{Participant identifier.}
  \item{architecture}{Processing architecture: \code{"parallel"} (the default and the existing behavior) or \code{"serial"}.}
  \item{serial_mode}{For a serial architecture, \code{"self-terminating"} stops as soon as the logical rule is determined, while \code{"exhaustive"} completes both channels. Serial self-termination is implemented for explicitly requested OR and AND rules; when \code{logical_rules} is omitted, serial self-termination defaults to OR. ID is not supported for serial self-termination.}
  \item{serial_order}{Serial channel order. \code{"A"} and \code{"B"} always check the named channel first; \code{"random"} samples the first channel independently on each trial using \code{p_A_first}; \code{"fixed"} is an alias for fixed A-first order.}
  \item{p_A_first}{For random serial order, the proportion of trials that check channel A first. Must be in [0, 1]; the default is 0.5.}
  \item{serial_first, serial_stopping, serial_p_A_first, serial_type}{Optional aliases for selecting the fixed first channel, serial stopping mode, and random-order probability, respectively.}
  \item{params, nSubjects, subject_ids, combine}{Hierarchical specifications for \code{simulate_sft_group}.}
  \item{vno, t, A, b, t0, sv, interval}{Parameters for LBA matching helper.}
  \item{backend, dt, max_time, bridge, bridge_exp_cutoff, adaptive, adaptive_mode, adapt_factor, adapt_eps, p_tol, curv_tol, kappa, tau, rho}{Simulation, integration, and backend controls for OU/LBA. \code{rho} is an alias for \code{tau}. In serial architecture, LBA capacity modifiers (including \code{kappa}, \code{tau}/\code{rho}, and capacity terms) and OU cross-channel interaction terms are set to zero because channels do not run concurrently.}
  \item{...}{Additional arguments passed to \code{simulate_sft}.}
}
\value{
  A comprehensive \code{data.frame} of simulated response times and accuracy, with columns 
  \code{RT}, \code{Channel1}, \code{Channel2}, \code{Correct}, \code{Subject}, and \code{Condition}. Serial outputs additionally record \code{Architecture}, \code{SerialMode}, \code{SerialOrder}, and \code{p_A_first}; \code{SerialOrder} is the channel checked first on serial trials. Parallel output columns remain unchanged.
}
\details{
  \code{simulate_sft_group} acts as a wrapper calling \code{simulate_sft} per subject. 
  It handles a subject-by-parameter matrix or a common parameter vector recycled across \code{nSubjects}. 
  The returned data frames can be immediately passed to standard SFT functions like 
  \code{capacityGroup} and \code{sicGroupBF}.

  Parallel processing remains the default. In serial self-terminating OR trials,
  the first target ends the trial and two absences require both channels. In serial
  self-terminating AND trials, the first absence ends the trial and two targets
  require both channels. Serial exhaustive trials always add the two channel
  processing times in the selected order.
}
