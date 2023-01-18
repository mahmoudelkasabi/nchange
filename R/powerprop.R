## proportions *******************************************************************************
# Helper function: power of detecting net change achieved using n
powerprop <- function (P1,P2, PXY, n, gamma, deff, alt, sig.level = sig.level)
{
  if (alt == "one.sided")
    za <- qnorm(1 - sig.level)
  if (alt == "two.sided")
    za <- qnorm(1 - sig.level/2)

  zb <- 1-pnorm(za - (sqrt(n*(P1-P2)^2)/(sqrt(deff*((P1 * (1-P1)) + (P2 * (1-P2)) - 2*gamma*(PXY - P1 * P2))))))
  list(Power= round(zb,2), n = ceiling(n), gamma = gamma)
}
