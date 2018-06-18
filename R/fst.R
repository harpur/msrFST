#' Calculate fst between two populations in an ms output
#'
#' This function allows you to load PFAM results to datafra,e
#' @param x the vector of alleles from ms
#' @param s1 pop size of the first population (first to appear in ms output)
#' @param s2 pop size of the second population (second to appear in ms output)
#' @keywords fst
#' @export
#' @examples
#' fst(x)


fst <- function(x, s1 = 24, s2 = 60){
  #http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0135368
  r <- 2 #always
  pop.size <- s1 + s2
  #allele freqs
  #pop.1.f <- freq(pop.1.sites)
  #pop.2.f <- freq(pop.2.sites)
  pop.1.f <- colSums(x[c(1:s1),]) / s1
  pop.2.f <- colSums(x[c((s1+1):pop.size),]) / s2
  p.bar <- (s1 * pop.1.f / pop.size) + (s2 * pop.2.f / pop.size)
  n.bar <- (s1 / r) + (s2 / r)
  #calculate S2
  S2 <- 1 / ((r - 1) * n.bar)
  temp <- (s1 * ((pop.1.f - p.bar)^2))
  temp2 <- (s2 * ((pop.2.f - p.bar)^2))
  S2 <- S2 * (temp + temp2)
  #calculate T1
  temp <- 1 / ((2 * n.bar) - 1 )
  temp1 <- (p.bar * (1 - (p.bar))) - (((r - 1)/r) * S2)
  T1 <- S2 - (temp * temp1)
  #calculate nc
  temp <- (1 / (r - 1))
  temp1 <- pop.size - ((s1^2 + s2^2) / (s1 + s2))
  nc <- temp * temp1
  #calculate T2
  temp <- ((2 * nc) - 1) / ((2 * n.bar) - 1)
  temp1 <- p.bar * (1 - p.bar)
  temp2 <- ((2 * (r - 1)) * (n.bar - nc) ) / ((2 * n.bar) - 1)
  temp3 <- S2 / r
  T2 <- (temp * temp1) + ((1 + temp2) * temp3)
  #output Fst
  fst.out <- T1/T2
  return(as.numeric(mean(fst.out)))

}
