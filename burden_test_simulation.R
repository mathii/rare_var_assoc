## Simulate a simple example of a burden test
## This simulates n independent loci and indendent traits, and
## for each trait simulates n.m snps. Then runs a burden test for
## the traits as in Morris & Zeggini. I.e. we are looking at n.loci
## different genes and traits. 

source("kpop.R")
source("coalescent.r")
source("association.R")
library("RColorBrewer")
dyn.load("kpop.so")

##################################################################################################################################
## Parameters


k <- 20                               #There are k^2 populations on a k*k lattice
n <- 2                                  #Sample this many individuals from each square. 
M <- 0.01                            #Migration rate

mu <- 0                                 #Mean, no reason not to be zero for quantitative trait. 
beta <- 0                               #Genotype effect

n.loci <- 10000                      #Simulate this many independant loci ("genes")
n.mutations <- c(1,3,10)                   #Simulate this many mutations ("snps") per locus
n.repeats <- 10                              #repeat simulations this many times and average the qq plots. 

min.f <- 0                              #Min and max frequency of variant to consider. 
max.f <- 0.04

################################################################################################################################## 

risk.funs <- c( hi.square.risk.fun, hi.gauss.blob.risk.fun )
risk.names <- c( "hi_square_risk_fun", "hi_gauss_risk_fun" )

################################################################################################################################## 

n.sim <- rep(n, k*k)

for( r in 1:length(risk.funs) ){
  risk.fun <- risk.funs[[r]]
  risk.name <- risk.names[r]
  summary.list <- list()
  cat(paste(Sys.time(), "- starting risk function", risk.name, "\n"))

  for(rp in 1:n.repeats){
    result.list <- list()
    cat(paste(Sys.time(), "- starting repeat", rp, "\n"))
    
    for(n.muts in n.mutations){
      cat( paste(Sys.time(), "- considering", n.muts, "rare variants per gene\n" ) )
      
      collapsed.genotypes <- matrix(0,nrow=n*k*k, ncol=n.loci)
      
      for( g in 1:n.loci ){
        ## Test one gene
        if(!(g %% 1000)){
          cat(paste(Sys.time(), "- starting locus", g, "\n"))
        }
      
        genotypes <- simulate.rare.genotypes(1, n.muts, n, k, M, min.f, max.f)
        collapsed.genotypes[,g] <- collapse.variants(genotypes)
      }
      
      cat(paste(Sys.time(), "- testing\n"))
    
      env <- get.env.risk(risk.fun, k, n.sim)
      ## a new quant trait for every locus. 
      quant.trait <- simulate.quant.trait(mu, collapsed.genotypes, beta, env=env, 0)
      result.list[[as.character(n.muts)]] <- lm.simple.vector(t(quant.trait), t(collapsed.genotypes))[,1]
    }
    
    summary.list[[rp]] <- plot.pvals.by.col(result.list,  minus.log10=TRUE, return.summary=TRUE, plot=FALSE)
  }

  summary <- Reduce("+",summary.list, 0)/length(summary.list)
  
  pdf(paste("burden_test_1210", risk.names[[r]], "_n_", n.loci , "_rp_", rp,".pdf", sep=""))
  par(lwd=2)
  plot.qq.by.af.from.summary(summary, feather = 4, feather.cex=0.4, feather.pch=20, legend.title=NULL)
  dev.off()

  write.table(summary, paste("burden test_1210", risk.name, "_n_", n.loci, "_rp_", n.repeats, ".txt", sep=""), row.names=FALSE, quote=FALSE)
}

    
