# powermeans: a helper function to find the power of detecting net change in means achieved using n
# Mahmoud Elkasabi

powermeans <- function (S2x, S2y, n, rho, gamma, deff, alt, del, sig.level = sig.level){
  if (alt == "one.sided")
    za <- stats::qnorm(1 - sig.level)
  if (alt == "two.sided")
    za <- stats::qnorm(1 - sig.level/2)

  zb <- 1-stats::pnorm(za - (sqrt(n*del^2)/(sqrt(deff*(S2x + S2y - 2*gamma*rho*sqrt(S2x)*sqrt(S2y))))))
  list(Power= round(zb,2), n = ceiling(n), gamma = gamma)
}
