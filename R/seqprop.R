## proportions *******************************************************************************
# Main function: Sequential assessment for power of detecting net change achieved using n
seqprop <- function (theta = 0.5, deff,
                     P1,P2, PXY, alt, sig.level = 0.05, power = 0.80,
                     P1.gross,P2.gross,PXY.gross,alt.gross, sig.level.gross = 0.05, pow.gross = 0.80)
{

  if (!(theta >= 0 & theta < 1))
    stop("theta must be in (0,1].\n")
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

  n11 <- nprop(P1.gross=P1.gross,P2.gross=P2.gross,PXY.gross=PXY.gross, deff=deff, alt.gross=alt.gross,
               sig.level.gross = sig.level.gross, pow.gross = pow.gross)$n1
  j = ((theta*n11) + theta - 1)/(1-theta)
  i <- 0

  n = ceiling(n11 + i + j) + 1

  power.check <- powerprop(P1=P1, P2=P2, PXY=PXY, n=n, gamma=n11/n, deff=deff, alt=alt, sig.level = sig.level)

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
      power.check <- powerprop(P1=P1, P2=P2, PXY=PXY, n=n, gamma=n11/n, deff=deff, alt=alt, sig.level = sig.level)

      samp1 <- rbind(samp1, round(power.check$n,0))
      samp2 <- rbind(samp2, round(power.check$Power,2))
      samp3 <- rbind(samp3, round(power.check$gamma,2))

    }

    RESULTS <- cbind.data.frame(samp1, samp2,samp3, n11, pow.gross,  row.names = NULL)

    names(RESULTS) <- c("n", "Power", "Gamma", "n11", "Power")
    RESULTS_last <- tail(RESULTS, n=1)
    names(RESULTS_last) <- c("n", "Power", "Gamma", "n11", "Power")

    structure(list("sequential power analysis" = RESULTS,
                   "final results" = RESULTS_last))
  }
}