#' Executes a sequential algorithm for sample size calculations for longitudinal surveys
#'
#' \code{seqmeans} executes a sequential algorithm for sample size calculation to measure net and gross changes concurrently in longitudinal surveys.
#' \code{seqmeans} deals with changes in means.
#'
#' @author Mahmoud Elkasabi.
#'
#' @param theta - 0 ≤ θ < 1 to set the lower bound of non-panel sample
#'
#' @param rho - unit-level correlation between x and y
#'
#' @param deff - design effect to adjust sample size for complex surveys
#'
#' @param S2x - unit variance of analysis variable x in sample t
#'
#' @param S2y - unit variance of analysis variable y in sample t+1
#'
#' @param alt - should the test of the net change be "1-sided" or "2-sided"
#'
#' @param del - size of the net change between the means to be detected
#'
#' @param sig.level - significance level of the hypothesis test of the net change
#'
#' @param power - desired power of the test of the net change
#'
#' @param S2o - common unit variance of analysis variable in panel sample
#'
#' @param alt.gross - should the test of the gross change be "1-sided" or "2-sided"
#'
#' @param del.gross - size of the gross change between the means to be detected
#'
#' @param sig.level.gross - significance level of the hypothesis test of the gross change
#'
#' @param pow.gross - desired power of the test of the gross change
#'
#' @import stats
#' @import utils
#'
#' @examples
#' # Calculate sample size of panel and fresh samples to estimate newt and gross changes
#'
#' sample_sequence <- seqmeans(
#' theta = 0.5,
#' rho = 0.9,
#' deff = 1,
#' S2x= 100,
#' S2y= 75,
#' alt= "one.sided",
#' del=3,
#' sig.level = 0.05,
#' power = 0.80,
#' S2o=200,
#' alt.gross="one.sided",
#' del.gross=3,
#' sig.level.gross = 0.05,
#' pow.gross = 0.80)
#'
#' @return panel sample size and overall sample size
#'
#' @export
seqmeans <- function (theta = 0.5, rho, deff,
                       S2x, S2y, alt, del, sig.level = 0.05, power = 0.80,
                       S2o, alt.gross, del.gross, sig.level.gross = 0.05, pow.gross = 0.80){

  if (!(theta >= 0 & theta < 1))
    stop("theta must be in (0,1].\n")
  if (rho < 0 | rho > 1)
    stop("rho must be in [0,1].\n")
  alt.ok <- alt %in% c("one.sided", "two.sided")
  if (!alt.ok)
    stop("alt must be either 'one.sided' or 'two.sided'.\n ")
  if (power < 0 | power > 1)
    stop("power must be in [0,1].\n")

  alt.ok <- alt.gross %in% c("one.sided", "two.sided")
  if (!alt.ok)
    stop("alt.gross must be either 'one.sided' or 'two.sided'.\n ")
  if (pow.gross < 0 | pow.gross > 1)
    stop("pow.gross must be in [0,1].\n")

  samp1 <- NULL
  samp2 <- NULL
  samp3 <- NULL

  n11 <- nmeans(S2o=S2o, rho=rho, deff=deff, alt.gross=alt.gross, del.gross=del.gross,
                  sig.level.gross = sig.level.gross, pow.gross = pow.gross)$n1
  j = ((theta*n11) + theta - 1)/(1-theta)
  i <- 0

  n = ceiling(n11 + i + j) + 1

  power.check <- powermeans(S2x=S2x, S2y=S2y, n=n, rho=rho, gamma=n11/n, deff=deff, alt=alt, del=del, sig.level = sig.level)

  if (power.check$Power >= power) {

    samp1 <- rbind(samp1, round(power.check$n,0))
    samp2 <- rbind(samp2, round(power.check$Power,2))
    samp3 <- rbind(samp3, round(power.check$gamma,2))

    RESULTS <- cbind.data.frame(samp1, samp2,samp3, n11, pow.gross,  row.names = NULL)
    names(RESULTS) <- c("n", "Power", "Gamma", "n11", "Power")

    structure(list("final results" = RESULTS))

  } else {

  while (power.check$Power < power){
    n = n+1
    i = i+1
    power.check <- powermeans(S2x=S2x, S2y=S2y, n=n, rho=rho, gamma=n11/n, deff=deff, alt=alt, del=del, sig.level = sig.level)

    samp1 <- rbind(samp1, round(power.check$n,0))
    samp2 <- rbind(samp2, round(power.check$Power,2))
    samp3 <- rbind(samp3, round(power.check$gamma,2))

  }

  RESULTS <- cbind.data.frame(samp1, samp2,samp3, n11, pow.gross,  row.names = NULL)

  names(RESULTS) <- c("n", "Power", "Gamma", "n11", "Power")
  RESULTS_last <- utils::tail(RESULTS, n=1)
  names(RESULTS_last) <- c("n", "Power", "Gamma", "n11", "Power")

  structure(list("sequential power analysis" = RESULTS,
                 "final results" = RESULTS_last))
  }
}
