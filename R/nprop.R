## proportions *******************************************************************************
# Helper function: Find n11 to measure gross change
nprop <- function (P1.gross=P1.gross,P2.gross=P2.gross,PXY.gross=PXY.gross,
                   deff, alt.gross, sig.level.gross, pow.gross)
{

  if (alt.gross == "one.sided")
    za <- qnorm(1 - sig.level.gross)
  if (alt.gross == "two.sided")
    za <- qnorm(1 - sig.level.gross/2)

  zb <- qnorm(1 - pow.gross)
  n1 <- deff * (((P1.gross * (1-P1.gross)) + (P2.gross * (1-P2.gross)) - 2 * (PXY.gross - P1.gross * P2.gross)) / ((P1.gross-P2.gross)^2)) * (za - zb)^2
  list(n1 = ceiling(n1))
}
