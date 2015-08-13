## Run a bunch of simulations for corrections, and then average the qq plots over them

source("kpop.R")
source("coalescent.r")
source("association.R")
library("RColorBrewer")

##################################################################################################################################
## Parameters


k <- 20                               #There are k^2 populations on a k*k lattice
n <- 2                                  #Sample this many individuals from each square. 
M <- 0.1                            #Migration rate

mu <- 0                                 #Mean, no reason not to be zero for quantitative trait. 
beta <- 0                               #Genotype effect

n.loci <- 1000                      #Simulate this many independant loci
n.mutations <- 10                   #Simulate this many mutations per locus
n.traits <- 10                         #Composite over this many trait simulations

################################################################################################################################## 
## Code starts here

dyn.load("kpop.so")

risk.funs <- c( hi.gauss.blob.risk.fun, hi.square.risk.fun)
risk.names <- c( "hi_gauss_risk_fun", "hi_square_risk_fun" )

################################################################################################################################## 

summaries<-list()

for( r in 1:length(risk.funs) ){
  risk.fun <- risk.funs[[r]]

  for(g in 1:n.traits){

    try({                          #Sometimes something fails (eg. mm doesn't converge - just ignore and try again in that case)
      cat(paste(Sys.time(), "- starting trait", g, "\n"))
      
      n.sim <- rep(n, k*k)
      
      genotypes <- simulate.genotypes(n.loci, n.mutations, n, k, M)
      
      env <- get.env.risk(risk.fun,k,n.sim)
      ## One trait to rule them all. 
      quant.trait <- simulate.quant.trait(mu, genotypes[,1,drop=FALSE], beta, env=env, 0)
      results <- t(apply( genotypes, 2, quant.trait.association.test, quant.trait=quant.trait ))
      results <- matrix(results, ncol=2 )
      
      cat( paste(Sys.time(), "- starting genomic control\n" ) )
      results.gc <- cbind(gc.correct(results[,1]),results[,2])
##       results.gcs <- cbind(gc.correct.by.bucket(results[,1], colMeans(genotypes)), results[,2])
      
      if(any("mm"==list.files())){
        cat( paste(Sys.time(), "- starting mm\n" ) )
        ks <- kinship(genotypes)
        results.mm <- mm.mm(quant.trait, genotypes, kinship=ks)
      }
      else{
        cat( paste(Sys.time(), "- couldn't find mm so no mixed models\n" ) )
      }
      
      cat( paste(Sys.time(), "- starting pca\n" ) )

      pca <- prcomp(t(genotypes), scale.=TRUE)
      pcs <- pca$rotation
      rm(pca)
      results.pca110 <- matrix(t(apply( genotypes, 2, quant.trait.association.test, quant.trait=quant.trait, covariates=pcs[,1:10])),ncol=2)
      
      cat( paste(Sys.time(), "- starting rare pca\n" ) )
      maf<-colMeans(genotypes)
      maf<-ifelse(maf<0.5, maf, 1-maf)
      genotypes.rare <- genotypes[,maf<0.04]
      pca.rare <- prcomp(t(genotypes.rare), scale.=TRUE) 
      pcs.rare <- pca.rare$rotation
      rm(pca.rare)
      results.pca.rare110 <- matrix(t(apply( genotypes, 2, quant.trait.association.test, quant.trait=quant.trait, covariates=pcs.rare[,1:10])),ncol=2)

      if(any("mm"==list.files())){
        results.list <-  list("Uncorrected"=results[,1], "GC"=results.gc[,1], "PCA 1:10"=results.pca110[,1], "Rare PCA 1:10"=results.pca.rare110[,1], "MM"=results.mm[,1] )
        else{
        results.list <-  list("Uncorrected"=results[,1], "GC"=results.gc[,1], "PCA 1:10"=results.pca110[,1], "Rare PCA 1:10"=results.pca.rare110[,1] )
      }          
      summaries[[g]] <- plot.pvals.by.col(results.list, minus.log10=TRUE, return.summary=TRUE, plot=FALSE)    

    })
      
  }

  summary <- Reduce("+",summaries, 0)/length(summaries)

  pdf(paste("composite_correction_", risk.names[[r]], "_t_", n.traits, "_n_", n.loci*n.mutations,".pdf", sep=""))
  par(lwd=2)
  plot.qq.by.af.from.summary(summary, feather = 4, feather.cex=0.4, feather.pch=20, legend.title=NULL)
  dev.off()

  write.table(summary, paste("composite_correction_", risk.names[[r]], "_t_", n.traits, "_n_", n.loci*n.mutations, ".txt", sep=""))
  cat(paste(length(summaries), "simulations succeded\n"))
}
