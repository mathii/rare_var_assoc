## Script to simualte random populations using Gil's kpop and coalescent scripts,
## and then do association testing on the results, also able to add in different
## effects according to some distribution of environmental risks.
source("kpop.R")
source("coalescent.r")
source("association.R")

##################################################################################################################################
## Parameters


k <- 20                                #There are k^2 populations on a k*k lattice
n <- 2                                #Get this many samples from each population
M <- 0.01                              #Migration rate
Ms <- c(0.1,0.5,1)

n.simulations=1000                      #Simulate this many coalescents
n.tests=1000                            #Simulate this many tests from each coalescent simulation

mu=0                                  #Baseline risk
beta=0                                #Genotype effect

tag <- "main"                           #tag for saving results plot

save.results=TRUE
load.saved.results = FALSE
plot.relative.inflation = FALSE


################################################################################################################################## 
## Code starts here

n.sim <- rep(n, k*k)

funs <-  c( gauss.blob.risk.fun, square.risk.fun)
names <- c( "gauss.risk", "square.risk" )

## Choose these to look at different distributions of risk
## funs <-  c( no.risk.fun, NW.risk.fun, square.risk.fun, six.square.risk.fun )
## names <- c("no.risk", "NW.risk", "square.risk", "six.square" )

## funs <-  c( square.risk.fun, big.square.risk.fun, big.big.square.risk.fun, big.big.big.square.risk.fun )
## names <- c("square", "big", "bigger", "biggest")

## funs <-  c( mid.square.risk.fun, mid.mid.square.risk.fun, mid.mid.mid.square.risk.fun )
## names <- c("smooth", "smoother", "smoothest")

## funs <-  c( no.risk.fun )
## names <- c("no_risk")

## funs <-  c( square.risk.fun, square.risk.fun, square.risk.fun, square.risk.fun )
## names <- c("very_slow", "slow", "fast", "very_fast")
## Ms <- c(0.01, 0.1, 0.5, 1)


## Load c extentions
dyn.load("kpop.so")

ng <- 2 + plot.relative.inflation + !load.saved.results

png(paste("quant_results_", tag, "_beta", beta, "_M_", M, "_k_", k, ".png", sep="" ), height=600 *length(funs), width=600*ng )
par(mfrow = c(length(funs),ng))

M.const=TRUE
if(M=="variable"){
  M.const=FALSE
}

for( r in 1:length(funs)){
  if(!M.const){
    M <- Ms[r]
  }
  risk.fun <- funs[[r]]
  if(load.saved.results){
    summary.inflation <- read.table(  paste( "saved.inflation.beta",beta, "M", M, names[r], "txt", sep="." ), sep ="\t", header = TRUE )
    summary.beta <- read.table(  paste( "saved.beta.beta", beta, "M", M, names[r], "txt", sep="." ), sep ="\t", header = TRUE )
  }
  else{
    env=get.env.risk(risk.fun,k,n.sim)
    results <- matrix(nrow=n.simulations*n.tests,ncol=3)
    for(i in 1:n.simulations){
      if(!(i%%1000)){
        cat( paste(Sys.time(), "- simulation", i, "\n" ) )
      }
      t<-simulate.tree.kpop.lattice.c(k=k, n=n.sim, M=M, annotate.tree=TRUE)
      results[((i-1)*n.tests+1):(i*n.tests),] <- run.n.simulations.vector(t, beta=beta, env=env, mu=mu, n.sim=n.tests )
    }
    colnames(results)<-c("p.vals", "beta", "afs")
    summary.inflation <- summarise.inflation.by.af(results[,1],results[,3],20,c(1e-4,1e-3, 1e-2, 0.1), maf=TRUE, log.scale=TRUE)
    summary.beta <- summarise.effect.size.by.af(results[,2], results[,1], results[,3], 20, c(0.0001, 0.001, 0.01, 0.1), TRUE, TRUE )
  }

  gc.correct(results[,1])
  
  plot.mutation(get.env.risk(risk.fun,k,1),k,rep(1,k), main=paste("M =", M))
  
  if(!load.saved.results){
    qq.sum <- plot.qq.by.af(results[,1], results[,3],c(0,0.04, 0.1, 0.5) ,maf=TRUE,xlim=c(0,6), ylim=c(0,10), main = "QQ plot by MAF", return.summary=TRUE)
  }

  plot.inflation.by.af(summary.inflation, main="P Values by quartile", xlab = "log10 Minor Allele Frequency", ylab = "Inflation in -log10 P value", ylim=c(0,4))
  if(plot.relative.inflation){
    if(r==1){baseline.inflation <- summary.inflation}
    relative.inflation <- summary.inflation
    relative.inflation[,3:ncol(relative.inflation)] <- summary.inflation[,3:ncol(relative.inflation)]- baseline.inflation[,3:ncol(relative.inflation)]
    plot.inflation.by.af(relative.inflation, main="P Values by quartile", xlab = "log10 Minor Allele Frequency", ylab = "Relative Inflation in -log10 P value", ylim=c(0,3))
  }
    
  if(!load.saved.results & save.results){
    write.table(summary.inflation, paste( "saved.inflation.beta",beta, "M", M, names[r],k, "txt", sep="." ), sep = "\t", row.names=FALSE, quote=FALSE )
    write.table(summary.beta, paste( "saved.beta.beta", beta, "M", M, names[r],k, "txt", sep="." ), sep = "\t", row.names=FALSE, quote=FALSE )
    write.table(qq.sum, paste( "saved.qq.beta", beta, "M", M, names[r],k, "txt", sep="." ), sep = "\t", row.names=FALSE, quote=FALSE )    
  }
}
dev.off()
