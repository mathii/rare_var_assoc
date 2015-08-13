## Functions to calculate allele sharing measures

source("kpop.R")
library(RColorBrewer)

## For a sample of k*k*n individuals ordered as in kpop.R, calculate the distance from
## individual i to each of the others. 

sample.distance <- function(i, n, k){
  i.pos <- decode.m(ceiling(i/n), k)
  all.pos <- sapply( ceiling((1:(n*k*k))/n), decode.m, k)
  d <- colSums(abs(all.pos-i.pos)) #Manhattan distance
  return(d)
}

sample.distance.matrix <- function(n,k){
  nkk <- n*k*k
  dm <- matrix(0,nrow=nkk,ncol=nkk)
  for( i in 1:(nkk) ){
    dm[i,] <- sample.distance(i, n, k)
  }
  return(dm)
}

## For rare, low frequency, and common alleles, compute the excess that are shared at each distance from
## each point. 

allele.sharing.by.distance <- function(gt, n, k){

  if(NROW(gt)!=n*k*k){stop("Genotype matrix is the wrong size")}

  n.gt <- NROW(gt)
  n.lc <- NCOL(gt)
  
  max.dist <- 2*k
  results <- matrix(0, ncol=max.dist-1, nrow=n.lc)

  a.f <- colMeans(gt)
  a.b <- ifelse(a.f < 0.04, 1, ifelse(a.f < 0.1, 2, 3))
  d.m <- sample.distance.matrix(n,k)+1  #add one for indexing
  
  for(m in 1:n.lc){
    this.shared <- rep(0,max.dist-1)    #How many alleles shared at each distance
    this.all <- rep(0,max.dist-1)       #How many alleles in total at each distance

    for(i in 1:n.gt){
      if(gt[i,m]==1 & i<n.gt){
        for(j in (i+1):n.gt){
          this.all[d.m[i,j]]=this.all[d.m[i,j]]+1
          if(gt[j,m]==1){
            this.shared[d.m[i,j]]=this.shared[d.m[i,j]]+1
          }
        }
      }
    }

    results[m,] <- this.shared/this.all
  }
  results <- results/a.f

  summary.results <- matrix(0,ncol=max.dist-1,nrow=3)
  for(i in 1:3){
    summary.results[i,] <- colMeans(results[a.b==c(1,2,3)[i],], na.rm=TRUE)
  }
  
  return(summary.results)
        
}

## Plot the excess allele sharing, assuming a three row matrix
## for rare, low frequency, and common alleles. 
plot.excess.allele.sharing <- function(sr, labels=c("<0.04", "0.04-0.1", ">0.1"), log.scale=FALSE, ... ){

  if(NROW(sr)!=length(labels)){stop("Please supply labels of the right length")}
  if(log.scale){sr <- log10(sr)}
  cols <- brewer.pal(NROW(sr), "Set1")

  yl="Excess allele sharing"
  if(log.scale){yl=expression(paste("lo",g[10],"(Excess allele sharing)",sep=""))}
  
  plot(1:length(sr[1,]), sr[1,], bty="n",type="l", col=cols[1], xlab="Distance", ylab=yl, ...)

  if(NROW(sr)>1){
    for(i in 2:NROW(sr)){
      lines(1:length(sr[1,]), sr[i,], col=cols[i])
    }            
  }

  abline(h=1-log.scale, lty=2, lwd=3)
  
  legend("topright", labels, col=cols, bty="n", lty=1, title="DAF")
  
}
