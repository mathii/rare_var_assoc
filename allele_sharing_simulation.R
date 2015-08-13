## Calculate allele sharing measures for our strubtured populations. 

source("kpop.R")
source("coalescent.r")
source("association.R")
source("allele_sharing_lib.R")
dyn.load("kpop.so")

##################################################################################################################################
## Parameters


k <- 20                                #There are k^2 populations on a k*k lattice
n <- 2                                 #Get this many samples from each population
M <- 10                              #Migration rate

n.loci=1000                      #Simulate this many coalescents
n.mutations=1                            #Simulate this many tests from each coalescent simulation

################################################################################################################################## 
## Script

n.sim <- rep(n, k*k)

gt <- simulate.genotypes(n.loci, n.mutations, n, k, M)

asd <- allele.sharing.by.distance(gt, n, k)

write.table(asd, paste("Allele_sharing_M_", M, "_k_", k, "_N_", n.loci*n.mutations, ".new.txt", sep=""), sep=",", col.names=FALSE, row.names=FALSE)

pdf(paste("Allele_sharing_M_", M, "_k_", k, "_N_", n.loci*n.mutations, ".pdf", sep=""))
par(lwd=4)
par(cex=1.4)
par(mar=c(4.2,4,2,2))
plot.excess.allele.sharing(asd, log.scale=TRUE, ylim=c(-2,2), xlim=c(0,2*k) )
dev.off()

