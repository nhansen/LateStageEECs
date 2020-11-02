setwd("/Users/nhansen/Bell_EEC")

# This analysis estimates the power to attain significance (q<0.1 for
# Benjamini-Hochberg FDRs) for p-values calculated by MutSigCV from
# 270 early stage tumors for 12 genes that had been discovered to be 
# significantly mutated in late-stage tumors. The calculation uses
# the binomial power model described in "Discovery and saturation
# analysis of cancer genes across 21 tumor types", Nature. 2014 Jan 23;
# 505(7484):495-501.
#
# Nancy F. Hansen
# October 30, 2020

# number of genes to be tested
tests <- 12
# number of early-stage tumors sequenced
N <- 270
# average gene length
L <- 1500
# gene-specific mutation rate valid for ~90% of genes
fg <- 3.9
# fraction of mutated tumors missed
m <- 0.1

# probability a tumor will have a non-synonymous mutation due to background
p_0 <- function(mu) {
  p <- 1 - (1 - mu*fg)^(3*L/4)
  return(p)
}

# probability a tumor will be mutated due to a gene's SMG status
p_1 <- function(r, mu) {
    p <- p_0(mu) + r * (1 - m)
    return(p)
}

# number of tumors that must be mutated to attain significance
Nmutmin <- function(mu) {
  pbg <- p_0(mu)
  perc <- 1 - 0.1/tests
  Nmut <- qbinom(perc, N, pbg)
  return(Nmut)
}

# probability that tumors mutated with p_1 will attain significance
binompower <- function(r, mu) {
  nminmut <- Nmutmin(mu)
  ptum <- p_1(r, mu)
  power <- pbinom(nminmut, N, ptum, lower.tail=FALSE)
  return(power)
}

# code to generate plot in Figure S##:
rvalues <- c(0.05, 0.06, 0.07, 0.08, 0.09, 0.1) # all plotted potential percentages of tumors with nonsyn mutations
allrates <- seq(1,20,0.01) # rate range for plotting (in 1 mutation/Mb)

minratevalues <- list() # will hold values of mut rate where power curve has minima
minpowercurve <- list() # will hold power values at minima

i <- 1
for (r in rvalues) { # generate a curve of local minima power for each possible r value:
  allpowercurve <- sapply(allrates, function(x) {binompower(r, x/1000000)})
  minpointindex <- 1
  minpowercurve[[i]] <- list()
  minratevalues[[i]] <- list()
  for (powerindex in seq(2, length(allpowercurve) - 1)) {
    if ((allpowercurve[powerindex-1] > allpowercurve[powerindex]) && (allpowercurve[powerindex+1] > allpowercurve[powerindex])) {
      minratevalues[[i]][[minpointindex]] <- allrates[powerindex]
      minpowercurve[[i]][[minpointindex]] <- allpowercurve[powerindex]
      minpointindex <- minpointindex + 1
    }
  }
  i <- i+1
}
curvecolors <- c("blue", "darkred", "brown", "orange", "purple", "darkgreen")
rvallabels <- lapply(rvalues, function(x) {paste("r=", x*100, "%", sep="")}) # labels for legend

plot(minratevalues[[1]], minpowercurve[[1]], type="l", ylab="Power to detect an SMG", xlab="Background mutation rate/Mb", col=curvecolors[[1]], ylim=c(0.8, 1.0), main="Power to detect early stage significantly mutated genes\n(270 tumors, 12 tests)")
points(minratevalues[[2]], minpowercurve[[2]], type="l", col=curvecolors[[2]])
points(minratevalues[[3]], minpowercurve[[3]], type="l", col=curvecolors[[3]])
points(minratevalues[[4]], minpowercurve[[4]], type="l", col=curvecolors[[4]])
points(minratevalues[[5]], minpowercurve[[5]], type="l", col=curvecolors[[5]])
points(minratevalues[[6]], minpowercurve[[6]], type="l", col=curvecolors[[6]])

legend(15, 0.88, lty=1, col=c("blue", "darkred", "brown", "orange", "purple", "darkgreen"), rvallabels)

