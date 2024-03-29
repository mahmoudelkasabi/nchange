# powerprop: a helper function to find the power of detecting net change in proportions achieved using n
# Mahmoud Elkasabi

powerprop <- function (P1,P2, PXY, n, gamma, deff, alt, sig.level = sig.level){
  if (alt == "one.sided")
    za <- stats::qnorm(1 - sig.level)
  if (alt == "two.sided")
    za <- stats::qnorm(1 - sig.level/2)

  zb <- 1-stats::pnorm(za - (sqrt(n*(P1-P2)^2)/(sqrt(deff*((P1 * (1-P1)) + (P2 * (1-P2)) - 2*gamma*(PXY - P1 * P2))))))
  list(Power= round(zb,2), n = ceiling(n), gamma = gamma)
}
