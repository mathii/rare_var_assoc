## Functions to simulate association testing for data generated with
## Gil's simulation scripts - kpop.R and coalescent.R
source("kpop.R")
## source("emma.R")

#########################################################
# environmental risk functions
#########################################################

no.risk.fun <- function(x,y){return(0*x)}
NW.risk.fun <- function(x,y){return(0.25*(y-x+1))}
N.risk.fun <- function(x,y){return(0.5*y)}
blob.risk.fun <- function(x,y){return(ifelse( sqrt((x-0.25)^2 + (y-0.75)^2) < 0.25, 0.5, 0))}
center.risk.fun <- function(x,y){return(ifelse( sqrt((x-0.5)^2 + (y-0.5)^2) < 0.25, 0.5, 0))}
big.square.risk.fun <- function(x,y){return(ifelse( y>0.1 & y <0.5 & x>0.5 & x < 0.9, 0.4, 0))}


square.risk.fun <- function(x,y){return(ifelse( y>0.2 & y <0.36 & x>0.6 & x < 0.76, 1, 0))}
hi.square.risk.fun <- function(x,y){return(ifelse( y>0.2 & y <0.36 & x>0.6 & x < 0.76, 2, 0))}
mid.square.risk.fun <- function(x,y){return(ifelse( y>0.2 & y <0.3 & x>0.6 & x < 0.7, 1, ifelse(y>0.15 & y < 0.35 & x > 0.54 & x < 0.76, 0.5, 0)))}
mid.mid.square.risk.fun<- function(x,y){return(ifelse( y>0.2 & y <0.25 & x>0.6 & x < 0.65, 1, ifelse(y>0.15 & y < 0.3 & x > 0.54 & x < 0.71, 0.7, ifelse(y>0.09 & y < 0.36 & x > 0.49 & x < 0.76, 0.35, 0))))}
mid.mid.mid.square.risk.fun<- function(x,y){return(ifelse( y>0.2 & y <0.25 & x>0.6 & x < 0.65, 0.8, ifelse(y>0.15 & y < 0.3 & x > 0.54 & x < 0.71, 0.6, ifelse(y>0.09 & y < 0.36 & x > 0.49 & x < 0.76, 0.4, ifelse(y>0.04 & y < 0.41 & x > 0.44 & x < 0.81, 0.2, 0)))))}
big.bad.square.risk.fun <- function(x,y){return(ifelse( y>0.2 & y <0.6 & x>0.35 & x < 0.76, 1, 0))}

big.square.risk.fun <- function(x,y){return(0.75*ifelse( y>0.2 & y <0.3 & x>0.6 & x < 0.7, 1, ifelse(y>0.15 & y < 0.35 & x > 0.54 & x < 0.76, 1, 0)))}
big.big.square.risk.fun<- function(x,y){return(0.6*ifelse( y>0.2 & y <0.25 & x>0.6 & x < 0.65, 1, ifelse(y>0.15 & y < 0.3 & x > 0.54 & x < 0.71, 1, ifelse(y>0.09 & y < 0.36 & x > 0.49 & x < 0.76, 1, 0))))}
big.big.big.square.risk.fun<- function(x,y){return(0.4285*ifelse( y>0.2 & y <0.25 & x>0.6 & x < 0.65, 1, ifelse(y>0.15 & y < 0.3 & x > 0.54 & x < 0.71, 1, ifelse(y>0.09 & y < 0.36 & x > 0.49 & x < 0.76, 1, ifelse(y>0.04 & y < 0.41 & x > 0.44 & x < 0.81, 1, 0)))))}

two.square.risk.fun <- function(x,y){return(ifelse((y>0.2 & y<0.36 & x>0.6 & x<0.76)|(y>0.7 & y<0.86 & x>0.1 & x<0.26) , 1, 0))}
three.square.risk.fun <- function(x,y){return(ifelse((y>0.2 & y<0.36 & x>0.6 & x<0.76)|(y>0.7 & y<0.86 & x>0.1 & x<0.26)|(y>0.2 & y<0.36 & x>0.1 & x<0.26) , 1, 0))}
four.square.risk.fun <- function(x,y){return(ifelse((y>0.2 & y<0.36 & x>0.65 & x<0.8)|(y>0.7 & y<0.86 & x>0.2 & x<0.36)|(y>0.2 & y<0.36 & x>0.2 & x<0.36)|(y>0.7 & y<0.86 & x>0.65 & x<0.8) , 1, 0))}

as.big.blob.risk.fun <- function(x,y){return(ifelse( sqrt((x-0.25)^2 + (y-0.75)^2) < 0.5, 0.5, 0))}
six.square.risk.fun <- function(x,y){return(ifelse( y>0.5 & x < 0.7, 0.5, 0))}
gauss.blob.risk.fun <- function(x,y){return(1*dnorm(sqrt((x-0.25)^2 + (y-0.75)^2)/0.25))}
hi.gauss.blob.risk.fun <- function(x,y){return(2*dnorm(sqrt((x-0.25)^2 + (y-0.75)^2)/0.25))}


#########################################################
# generate a bunch of random genotypes using the structured
## coalescent. We generate n.loci independent trees and
## on each tree genrate n.mutions independent mutations
## The population lives on a k*k grid with n individuals
## in each square. The migration rate is M. 
#########################################################

simulate.genotypes <- function(n.loci, n.mutations, n, k, M, n.sim=rep(n,k*k)){

  n.ind <- sum(n.sim)
  genotypes <- matrix(n.ind*n.loci*n.mutations, nrow=n.ind, ncol=n.loci*n.mutations, )
  
  for( i in 1:n.loci ){
    if(!(i%%1000)){
      cat( paste(Sys.time(), "- simulation", i, "\n" ) )
    }
    t<-simulate.tree.kpop.lattice.c(k=k, n=n.sim, M=M, annotate.tree=TRUE)
    mutations <-  choose.mutation(t$tree, TRUE, n.mutations=n.mutations)
    n.indiv <- (nrow(t$tree)+1)/2
    all.descendants <-find.all.descendants(t$tree)
    these.genotypes <- sapply(all.descendants[mutations], indices.to.genotype, n.indiv)
    genotypes[,((i-1)*n.mutations+1):(i*n.mutations)] <- these.genotypes
  }
  return(genotypes)
}

#########################################################
# generate a bunch of random genotypes using the structured
## coalescent. As simualte genotpyes, but keep only those
## with frequency in the range given. If n.sim is specified
## then n is ignored. 
#########################################################

simulate.rare.genotypes <-  function(n.loci, n.mutations, n, k, M, min.f=0, max.f=0.05, n.sim=rep(n, k*k)){
  n.ind <- sum(n.sim)
  genotypes <- matrix(n.ind*n.loci*n.mutations, nrow=n.ind, ncol=n.loci*n.mutations, )
  
  for( i in 1:n.loci ){

    if(!(i%%1000)){
      cat( paste(Sys.time(), "- simulation", i, "\n" ) )
    }

    this.locus.genotypes <- matrix(-1, nrow=n.ind, ncol=n.mutations, )
    n.found <- 0                        #Number of good mutations generated. 
    t<-simulate.tree.kpop.lattice.c(k=k, n=n.sim, M=M, annotate.tree=TRUE)

    ## Some malformed trees need investigating. 
    if(any(is.na(t$tree[,6]))){
      while(any(is.na(t$tree[,6]))){
        t<-simulate.tree.kpop.lattice.c(k=k, n=n.sim, M=M, annotate.tree=TRUE)
      }
    }
    
    while(n.found<n.mutations){
      mutations <-  choose.mutation(t$tree, TRUE, n.mutations=n.mutations)
      n.indiv <- (nrow(t$tree)+1)/2
      all.descendants <-find.all.descendants(t$tree)
      these.genotypes <- sapply(all.descendants[mutations], indices.to.genotype, n.indiv)
      MAF <- colMeans(these.genotypes)
      MAF <- ifelse(MAF<0.5, MAF, 1-MAF)
      in.range <- MAF>=min.f & MAF<=max.f
      if(sum(in.range>0)){
        n.required <- n.mutations-n.found
        n.used <- min(n.required,sum(in.range))
        this.locus.genotypes[,(n.found+1):(n.found+n.used)] <- these.genotypes[,in.range,drop=FALSE][,1:n.used,drop=FALSE]
        n.found <- n.found+n.used
      }
    }
    
    genotypes[,((i-1)*n.mutations+1):(i*n.mutations)] <- this.locus.genotypes
  }
  return(genotypes)
}

#########################################################
# Generate random case control status for samples,
# given (haploid) genotyes, genotype risk, environmental
# risk, in log odds units and GE interaction term.
#########################################################

simulate.cc.status <- function(mu, genotypes, beta=0, env=0, interact=0){
  log.odds <- mu + beta*genotypes + env + interact*genotypes*env
  status <- rbinom(genotypes, 1, exp(log.odds)/(1+exp(log.odds)) )
  return(status)
}

#########################################################
# Do an association test for case control data, and
# return a p.value. 
#########################################################

cc.association.test <- function(cc.status, genotypes){
  model <- glm(cc.status ~ genotypes, family="binomial")
  p.value <- summary(model)$coefficients["genotypes","Pr(>|z|)"]
  beta<-summary(model)$coefficients["genotypes","Estimate"]  
  return(c(p.value, beta))
}

#########################################################
# Generate quantitative trait values, given genotypes,
# additive genetic effects, environment etc.. Variance
# of trait value is 1
#########################################################

simulate.quant.trait <- function(mu, genotypes, beta=0, env=0, interact=0){
  mean = mu + beta*genotypes + env + interact*genotypes*env
  if(is.null(dim(genotypes))){
    trait <- rnorm(genotypes, mean, sd=1 )
  }
  else{
    trait <- matrix(rnorm(prod(dim(genotypes)), mean, sd=1 ),dim(genotypes)[1], dim(genotypes)[2])
  }
  return(trait)
}

#########################################################
# Do an association test for a quantitative trait, and
# return a p.value. use lm.simple instead, since that
# will be faster for the simple calculation we need here
# but doesn't allow for extra covariates
#########################################################

quant.trait.association.test <- function(quant.trait, genotypes, covariates=NULL){
  if(is.null(covariates)){
    model <- lm(quant.trait ~ genotypes)
  }
  else{
    model <- lm(quant.trait ~ genotypes + covariates)
  }
  p.value <- summary(model)$coefficients["genotypes","Pr(>|t|)"]
  beta<-summary(model)$coefficients["genotypes","Estimate"]
  return(c(p.value, beta))
}

#########################################################
# simple linear model
# return p value and beta for y~x, 1-D
## for lm.simple.vector you probably want to transpose
## the genotypes and traits, so that nrow=n.snps and
## ncol=n.individuals. 
#########################################################

lm.simple <- function(y,x){
  x.bar<-mean(x)
  x2.bar<-mean(x*x)
  y.bar<-mean(y)
  xy.bar<-mean(x*y)
  x.var<-x2.bar-x.bar*x.bar
  df=length(x)-2
  
  b.hat<-(xy.bar-x.bar*y.bar)/x.var
  a.hat<-y.bar-b.hat*x.bar

  e2.hat<-y-a.hat-b.hat*x

  sd.b.hat<-sqrt( (sum(e2.hat*e2.hat)) / (x.var*df*(df+2)) )

  p <- 2*pt(-abs(b.hat/sd.b.hat),df=df)

  return(c(p,b.hat))
}

lm.simple.vector <- function(y,x){
  x.bar<-rowMeans(x)
  x2.bar<-rowMeans(x*x)
  y.bar<-rowMeans(y)
  xy.bar<-rowMeans(x*y)
  x.var<-x2.bar-x.bar*x.bar
  df=dim(x)[2]-2
  
  b.hat<-(xy.bar-x.bar*y.bar)/x.var
  a.hat<-y.bar-b.hat*x.bar

  e2.hat<-y-a.hat-b.hat*x
  sd.b.hat<-sqrt( (rowSums(e2.hat*e2.hat)) / (x.var*df*(df+2)) )
  
  p <- 2*pt(-abs(b.hat/sd.b.hat),df=df)

  return(cbind(p,b.hat))
}

lm.simple.c <- function(y,x){
  n <- length(x)
  df<-n-2
  results<-c(0,0)
  res <- .C("lm_simple_c", as.double(y), as.double(x), as.integer(n), as.double(results))
  results <- res[[4]]
  results[1] <- 2*pt(-abs(results[2]/results[1]),df=df)
  return(results)
}

#########################################################
# Get risk. Given a function of x,y positions between
# 0 and 1, get the risk for the populations on the k*k
# lattice. 
#########################################################

get.env.risk <- function( risk.function, k, n.sim ){
  points <- (sapply(1:(k*k), decode.m, k=k )-1)/(k-1)
  risk <- risk.function(points[1,],points[2,] )
  env.risk <- rep(risk, times=n.sim)
  return( env.risk )
}

get.matrix.env.risk <- function(risk.function, k){
  return(matrix(mapply(risk.function, rep(((1:k)-0.5)/k, times=k), rep(((1:k)-0.5)/k, each=k ) ),nrow=k, ncol=k))
}

#########################################################
# Simulate risk, with some random genotypic effects.
# the effects are drawn from a normal
# distribution with mean and variance given by effect.params. 
#########################################################

get.genotypic.risk <- function(genotypes, effect.params=c(0,0.1)){
  n.snps <- NCOL(genotypes)
  snp.effects <- rnorm(n.snps, mean=effect.params[1], sd=effect.params[2])
  if(n.snps>1){
    risk <- colSums(apply(genotypes, 1, '*', snp.effects))
  }
  else{
    risk <- c(genotypes*snp.effects)
  }
  return(risk)
}

#########################################################
# Summarise the inflation in p-value by quantile for a
# list of p values, and allele frequencies. 
#########################################################

summarise.inflation.by.af <- function(p.vals, af, bins=20, quants = 0.5, maf=FALSE, log.scale=FALSE){

  if(log.scale){
    if(maf){range=c(-2+log10(0.5), log10(0.5))}
    else{range=c(-2,0)}
  }
  else{
    if(maf){range=c(0,0.5)}
    else{range=c(0,1)}
  }

  if(maf){af <- ifelse(af<0.5,af,1-af)}
  if(log.scale){af<-log10(af)}
  
  step = (range[2]-range[1])/bins
  bin.starts = range[1]+step*(0:(bins-1))
  bin.ends = bin.starts+step

  results<-data.frame(bin.starts=bin.starts, bin.ends=bin.ends)
  
  for(quant in quants){
    inf.factor <- c()
    for( b in 1:bins ){
      these.p.vals <- p.vals[(af>bin.starts[b])&(af<bin.ends[b])]
      quant.p <- quantile(these.p.vals,quant, na.rm=T)
      inf.factor <- c(inf.factor, -log10(quant.p)+log10(quant))
    }
    results[,as.character(quant)]<-inf.factor
  }

  return(results)
}

summarise.effect.size.by.af <- function(eff, p.vals, af, bins=20, quants=0.5, maf=TRUE, log.scale=FALSE ){
  if(log.scale){
    if(maf){range=c(-2+log10(0.5), log10(0.5))}
    else{range=c(-2,0)}
  }
  else{
    if(maf){range=c(0,0.5)}
    else{range=c(0,1)}
  }

  if(maf){af <- ifelse(af<0.5,af,1-af)}
  if(log.scale){af<-log10(af)}

  step = (range[2]-range[1])/bins
  bin.starts = range[1]+step*(0:(bins-1))
  bin.ends = bin.starts+step

  results<-data.frame(bin.starts=bin.starts, bin.ends=bin.ends)
  
  for(quant in quants){
    inf.factor <- c()
    for( b in 1:bins ){
      these.effects <- eff[(af>bin.starts[b])&(af<bin.ends[b])&(p.vals<quant)]
      inf.factor <- c(inf.factor, mean(these.effects))
    }
    results[,as.character(quant)]<-inf.factor
  }

  return(results)
}

summarise.statistic.by.af <- function(stat, af, lev, bins=20, levels=0.5, maf=TRUE, log.scale=FALSE){

  if(log.scale){
    if(maf){range=c(-2+log10(0.5), log10(0.5))}
    else{range=c(-2,0)}
  }
  else{
    if(maf){range=c(0,0.5)}
    else{range=c(0,1)}
  }

  if(maf){af <- ifelse(af<0.5,af,1-af)}
  if(log.scale){af<-log10(af)}

  step = (range[2]-range[1])/bins
  bin.starts = range[1]+step*(0:(bins-1))
  bin.ends = bin.starts+step

  results<-data.frame(bin.starts=bin.starts, bin.ends=bin.ends)
  
  for(level in levels){
    inf.factor <- c()
    for( b in 1:bins ){
      these.effects <- stat[(af>bin.starts[b])&(af<bin.ends[b])&(lev==level)]
      inf.factor <- c(inf.factor, mean(these.effects))
    }
    results[,as.character(level)]<-inf.factor
  }

  return(results)
}

#########################################################
# Combine inflation summaries - averaging quantiles
# is not quite the right thing to do, but better than
# nothing. 
#########################################################

combine.inflation.summaries <- function(sum1, sum2){
  if(!all(colnames(sum1)==colnames(sum2))){stop("Column names not equal")}
  if(!all(dim(sum1)==dim(sum2))){stop("Summary sizes not equal")}
  if(dim(sum1)[2]<3){stop("no data in summary")}
  if(!all(colnames(sum1)[1:2]==c("bin.starts", "bin.ends"))){stop("Wrong colnames")}
  if(!all(sum1[,1]==sum2[,1])|!all(sum1[,2]==sum2[,2])){stop("Different bins")}

  sum.new <- sum1
  n <- dim(sum1)[2]
  sum.new[,3:n] <- (sum1[,3:n]+sum2[,3:n])/2
  return(sum.new)
}
  
#########################################################
# Plot inflation by allele frequency
# given a summary from summarise.inflation.by.af
# allele frequencies, plot the inflation of
# p values at any particular quantile. 
#########################################################

plot.inflation.by.af <- function(inflation.summary, xlab="Allele Frequency", ylab=expression(paste("Inflation in lo",g[10]," P-value",sep="")), ylim=NA, legend.title = "Quantile", cols=get.cols(dim(inflation.summary)[2]-2), ...){

  range <- c(min(inflation.summary$bin.starts), max(inflation.summary$bin.ends))
             
  i=1
  quants<-colnames(inflation.summary)[3:ncol(inflation.summary)]
  bin.mids <- (inflation.summary$bin.starts+inflation.summary$bin.ends)/2
  for(quant in quants){
    inf.factor <-inflation.summary[,quant]

    if(quant==quants[1]){
      if(any(is.na(ylim))){
        ylim = c(0,max(inf.factor, na.rm=TRUE)*1.05)
      }
      xlim = c(floor(2*range[1])/2,ceiling(2*range[2])/2)
      
      plot( bin.mids, inf.factor,col = cols[i], xlim=xlim , ylim=ylim, xlab=xlab, ylab=ylab, type = "l", bty = "n", pch=20, cex=1/1.4,... )
    }
    else{
      points( bin.mids, inf.factor,col = cols[i], type = "l", pch=20, ... )
    }
    i=i+1
  }

  legend("topright", quants, col=cols, bty="n", lty=1, title=legend.title )
}


#########################################################
# QQ plots with a different line for different allele
# frequencies.
#########################################################

plot.qq.by.af <- function(p.vals, af, breaks=10, maf=FALSE, return.summary=FALSE, legend.title="MAF", ...){
  if(maf){
    af <- ifelse(af<0.5,af,1-af)
  }

  ## If we've just given a number then use that as the number of breaks
  if(length(breaks)==1){
    breaks=(0:breaks)/10/(1+maf)
  }
  rlist <- vector("list",length(breaks)-1)
  
  for( b in 1:(length(breaks)-1) ){
    rlist[[b]] <-  p.vals[af>breaks[b] & af<=breaks[b+1]]
    names(rlist)[b] <- paste( breaks[b], "-", breaks[b+1], sep = " ")
  }
  return(plot.pvals.by.col(rlist, minus.log10=TRUE, legend.title=legend.title, return.summary=return.summary, ... ))
}

#########################################################
## Plot qq from summary - plots the results of plot.qq.by.ag
## should look the same as the original
#########################################################

plot.qq.by.af.from.summary <- function(summary, legend.title="MAF",xlab = expression(paste("Expected -lo",g[10]," P",sep="")), ylab=expression(paste("Observed -lo",g[10]," P",sep="")), feather=Inf, feather.cex=1, feather.pch=20, cols=get.cols(length(names(summary))/2), ...){

  n.lines <- length(names(summary))/2

  for( i in 1:n.lines ){
    xpts <- summary[,(2*i-1)]
    these.p.vals <- summary[,2*i]

    xpts.line <- xpts[xpts<feather]
    these.p.vals.line <- these.p.vals[xpts<feather]

    xpts.dots <- xpts[xpts>=feather]
    yi <- these.p.vals[xpts>=feather]

    n <- length(yi)
    if(n>2){
      select <- (yi[1:(n-2)]+yi[3:n]-2*yi[2:(n-1)])>1e-10 #points where the line changes direction
      select <- c(FALSE, select, TRUE)
      xpts.dots <- xpts.dots[select]
      these.p.vals.dots <- yi[select]
    }
    if(i==1){
      plot(xpts.line, these.p.vals.line,col=cols[i],xlab = xlab, ylab= ylab, type = "l", bty = "n",...)
    }
    else{
      lines(xpts.line, these.p.vals.line, col=cols[i], type = "l", ...)
    }
    if(n>2){
      points(xpts.dots, these.p.vals.dots, col=cols[i], pch=feather.pch, cex=feather.cex, lwd=0.1)
    }
  }
  abline(0,1,col = "black", lty = 2)
  
  legend("bottomright", legend = names(summary[c(FALSE,TRUE)]), lty=1, col = cols, bty = "n", title=legend.title )
}


#########################################################
# QQ plots with a different line for columns. p.vals should
# be a list with one entry for each line that you
# want printed.
# set minus.log10 to TRUE if you want to plot -log10 p
# instead. 
#########################################################

plot.pvals.by.col <- function(p.vals, minus.log10=FALSE, type = "l", legend.title=NULL, plot=TRUE, return.summary=FALSE, ...){
  cols <- get.cols(length(p.vals))

  summary.grid <- 1000
  summary <- data.frame(matrix(0,nrow=summary.grid+1, ncol=2*length(p.vals)))
  
  for( i in 1:length(p.vals) ){
    these.p.vals <- sort(p.vals[[i]])
    N <- length(these.p.vals)
    xpts <- 1:N/(N+1)
    if(minus.log10){
      these.p.vals <- -log10(these.p.vals)
      xpts <- -log10(xpts)
    }
    if(i==1 & plot){
      plot(xpts, these.p.vals,col=cols[i],xlab = "Expected", ylab= "Observed", type = type, bty = "n",...)
    }
    else if(plot){
      lines(xpts, these.p.vals, col=cols[i], type = type, ...)
    }
    interpolate <- approx(xpts, these.p.vals, n=summary.grid+1)
    summary[,2*i-1] <- interpolate$x
    summary[,2*i] <- interpolate$y
    colnames(summary)[(2*i-1):(2*i)] <- c(paste(names(p.vals)[i], ".x", sep=""), names(p.vals)[i])
  }
  if(plot){
    abline(0,1,col = "black", lty = 3)
    legend("bottomright", legend = names(p.vals), lty=1, col = cols, bty = "n", title=legend.title )
  }
  if(return.summary){return(summary)}
  else{return()}
}

#########################################################
# Given an annotated tree, run a case/control simulations
# from it, with uniform mutations 
#########################################################

run.one.simulation <- function(t, beta=0, env=0, mu=0, case.control=FALSE ){
  genotypes <- add.mutation(t$tree, choose.mutation(t$tree, FALSE))
  if(case.control){
    cc.status <- simulate.cc.status(mu, genotypes, beta, env=env, 0)
    results <- cc.association.test(cc.status, genotypes )
  }
  else{
    quant.trait<- simulate.quant.trait(mu, genotypes, beta, env=env, 0)
    results <- quant.trait.association.test(quant.trait, genotypes)
  }
  afs <- sum(genotypes)/length(genotypes)
  return(list(p.vals=results[1], beta=results[2] ,afs=afs))
}

#########################################################
# Get genotypes, assuming we have first calculated a list
# of all the daughters for the tree. 
#########################################################

get.genotypes <- function(node, all.descendants ){
  	n<-(length(all.descendants)+1)/2;
  	seq<-rep(0, n)
	seq[all.descendants[[node]]]<-1
	return(seq)
}

#########################################################
# Given an annotated tree, run n case/control simulations
# from it, with uniform mutations. If genetic effects=TRUE
# then we simulate random genetic effects using loci
# independent to the simulated genotypes. genetic params
# are (n.loci, n, K, M mean, sd, min.f, max.f)
#########################################################

run.n.simulations <- function(t, beta=0, env=0, mu=0, case.control=FALSE, n.sims=1, genetic.risk=FALSE, genetic.params=c(100,2,20,0.01,0,0.5,0,0.1) ){
  mutations <-  choose.mutation(t$tree, TRUE, n.mutations=n.sims)

  overall.results<-matrix(nrow=n.sims, ncol=3)
  all.descendants <-find.all.descendants(t$tree)
  
  for(i in 1:n.sims){
    genotypes <- get.genotypes(mutations[i], all.descendants )
    if(genetic.risk){
      effect.gts <- simulate.rare.genotypes(genetic.params[1], 1, genetic.params[2], genetic.params[3], genetic.params[4], genetic.params[5], genetic.params[6])
      risk.gts <- get.genotypic.risk(effect.gts, genetic.params[7:8])
      env <- env+risk.gts
    }
      
    if(case.control){
      cc.status <- simulate.cc.status(mu, genotypes, beta, env=env, 0)
      results <- cc.association.test(cc.status, genotypes )
    }
    else{
      quant.trait<- simulate.quant.trait(mu, genotypes, beta, env=env, 0)
      results <- lm.simple(quant.trait,genotypes)
    }
    overall.results[i,] <- c(results, sum(genotypes)/length(genotypes) )
  }
  
##   overall.results<-data.frame(overall.results)
  names(overall.results) <- c("p.vals", "beta", "afs")
  return(overall.results)
}

indices.to.genotype <- function(indices, n){
  gt<-rep(0,n)
  gt[indices] <- 1
  return(gt)
}

run.n.simulations.vector <- function(t, beta=0, env=0, mu=0, case.control=FALSE, n.sims=1, genetic.risk=FALSE, genetic.params=c(100,2,20,0.01,0,0.5,0,0.1) ){  
  mutations <-  choose.mutation(t$tree, TRUE, n.mutations=n.sims)
  n.indiv <- (nrow(t$tree)+1)/2
  
  overall.results<-matrix(nrow=n.sims, ncol=3)
  all.descendants <-find.all.descendants(t$tree)
  genotypes <- sapply(all.descendants[mutations], indices.to.genotype, n.indiv)

  if(genetic.risk){
    effect.gts <- simulate.rare.genotypes(genetic.params[1], 1, genetic.params[2], genetic.params[3], genetic.params[4], genetic.params[5], genetic.params[6])
    risk.gts <- get.genotypic.risk(effect.gts, genetic.params[7:8])
    env <- env+risk.gts
  }
  
  if(case.control){
    stop("case control simualations not implemented in vectorised form")
  }
  else{
    quant.trait <- simulate.quant.trait(mu, genotypes, beta, env=env, 0)
##     results <- sapply(1:ncol(genotypes), function(i){ lm.simple(quant.trait[,i],genotypes[,i])})
    results <- lm.simple.vector(t(quant.trait), t(genotypes))
    results <- cbind(results, colMeans(genotypes))
  }
    
  results<-results
  return(results)
}

#########################################################
# Make a pretty qq-plot
#########################################################

qq.plot.of.p.values <- function( pvals, gi = NULL, add=FALSE, BW=FALSE, cex.factor=1, linecol="red", linewd=1, linety=2, ... ){
  p <- sort(unique(pvals))
  e = -log10((1:length(p))/length( p ))
  o = -log10(p)
  if(add){
    points( e, o, pch= 20, cex=0.5*cex.factor,... )
  }

  else{

    lo <- -log10(sapply(1:length(p),function(x) qbeta(.025,x,length(p)-x+1)))
    up <- -log10(sapply(1:length(p),function(x) qbeta(.975,x,length(p)-x+1)))

    
    plot( e, o, xlab = "Expected -log10(p-values)", ylab = "Observed -log10(p-values)", pch= 20, bty = "n", cex=0.3*cex.factor, cex.lab=cex.factor, cex.axis=cex.factor, xaxt="n", yaxt="n", type="n",... )
    axis(1, lwd=par("lwd"))
    axis(2, lwd=par("lwd"))
    if(!BW){
      polygon(c(e,rev(e)),c(lo,rev(up)),col="grey95",border="grey95",)
      lines(e, lo, col="grey")
      lines(e, up, col="grey")
    }
    points( e, o, pch= 20, cex=0.5*cex.factor, ... )
  }

  if(BW){
    lines(e,lo,lty=2)
    lines(e,up,lty=2)
    abline(0,1,col="black", lty=2)
  }
  else{
    abline(0,1,col=linecol, lwd=linewd, lty=linety)
  }
  
  if( !is.null(gi)){
    medx = median(e)
    medy = median(o)
    offset = 0.125 * max(e)
    lines(c(medx,medx-0.5*offset),c(medy,medy+offset),col="red",lty=3)
    text(medx-0.5*offset,medy+offset+0.1,substitute(lambda == gi, list(gi=round(gi,2)) ) )
  }
}

#########################################################
## Calculate Moran's I statistic to measure the
## spatial autocorrelation of the distirbution
## of a variant. Genotypes is a k*k vector and weights
## is a (k*k, k*k) matrix which tells us whether two
## nodes are adjacent. 
#########################################################

moran.I <- function(k, genotypes, weights){
  x <- genotypes-mean(genotypes)
  cov <- sapply(x, "*", x )
  a <- sum(weights*cov)
  var <- sum(x*x)
  I <- (k*k)*a/(sum(weights)*var)
  return(I)
}

#########################################################
## calculate Fst. Takes a genotype matrix with nrow
## individuals and ncol loci and a vector mapping
## individuals to 1:m populations. Haploid genotypes?
#########################################################

Fst <- function(genotypes, populations){
  if(!all(sort(unique(populations))==1:max(populations))){
    stop("Populations must range from 1:m, where m is the number of populations")
  }
  else{
    npops <- max(populations)
  }

  PI.within <- c()
  
  for(p in 1:npops){
    PI.within <- c(PI.within, sample.PI(genotypes[populations==p,]))
  }
  
  PI.within <- mean(PI.within)
  PI.between <- sample.PI(genotypes)

  fst <- 1-PI.within/PI.between
  return(fst)
}

#########################################################
## average number of pairwise differences between two
## samples. (or within one sample if only 1 provided)
#########################################################

sample.PI <- function(genotypes1, genotypes2=NULL){
  PI <- 0

  n <- nrow(genotypes1)
  if(is.null(genotypes2)){
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        PI <- PI+sum(abs(genotypes1[i,]-genotypes1[j,]))
      }
    }
    PI <- 2*PI/(n*(n-1))
  }
  else{
    m <- nrow(genotypes2)
    for(i in 1:n){
      for(j in 1:m){
        PI <- PI+sum(abs(genotypes1[i,]-genotypes2[j,]))
      }
    }
    PI <- PI/(m*n)
  }
  return(PI)
}
  
#########################################################
## Get n colors for coloring levels on the trees
#########################################################
get.cols <- function(n){
  if(!require(RColorBrewer)){
    return(rainbow(n))
  }
  else{
    if(n<9){return(brewer.pal(9,"Set1")[1:n])}
    else{return(rainbow(n))}
  }
}

#########################################################
## GC correction for p values
#########################################################

gc.correct <- function(p.vals){
  chi<- qchisq(p.vals, df=1, lower.tail=FALSE)
  gi <- median(chi)/qchisq(0.5, df=1, lower.tail=FALSE)
  if(gi < 0 ){
    cat("gi<0, not correcting\n")
  }
  else{
    cat(paste("gi=", gi, "\n", sep=""))
    chi <- chi/gi
    p.vals <- pchisq(chi, df=1, lower.tail=FALSE)
  }
  return(p.vals)
}

#########################################################
## GC correction by allele frequency buckets
#########################################################

gc.correct.by.bucket <- function(p.vals, freqs, buckets=c(0, 0.04, 0.1, 1)){
  for(i in 2:length(buckets)){
    include <- freqs>buckets[i-1] & freqs<buckets[i]
    p.vals[include] <- gc.correct(p.vals[include])
  }
  return(p.vals)
}

#########################################################
## Plot the positions of the individuals and color
## according to the values of a principal component
#########################################################

plot.pca <- function(pc, real.x, real.y, col.scheme="RdBu", ...){

  cols <- brewer.pal(8, col.scheme)
  pc <- (pc-min(pc))/(max(pc)-min(pc))
  pc <- floor(1 + 7*pc)

  plot(jitter(real.x), jitter(real.y), col=cols[pc], pch=20, ylab="N-S", xlab = "E-W", yaxt="n", xaxt="n", bty="o", ... )

}  

#########################################################
## An example plot showing generation of qt's
#########################################################

qt.example <- function(){
  plot(dnorm,xlim=c(-4,6), bty="n", xlab="Mean value", ylab="density", lty=2)
  cols <- brewer.pal(9, "YlOrRd")
  x <- seq(-4,6,0.01)
  for(i in 1:9){
    f<-function(x){dnorm(x-i/9)}
    y <- f(x)
    lines(x, y, col=cols[i])
  } 
}
#########################################################
## haploid kinship estimation 
#########################################################

haploid.kinship <- function(genotypes){
  N <- nrow(genotypes)
  L <- ncol(genotypes)
  K <- matrix(0,ncol=N,nrow=N)

  for(l in 1:ncol(genotypes)){
    K=K+haploid.kinship.1col(genotypes[,l,drop=FALSE])
  }

  K <- K/N
  return(K)
}

haploid.kinship.1col <- function(col){
  mc <- mean(col)
  k <- (col-mc)%*%t(col-mc)/(mc*(1-mc))
  return(k)
}

#########################################################
##  estimate relatedness matrix from matrix of who your
## nearest neighbout is at each locus. 
#########################################################

nn.relatedness <- function(nn){
  n <- dim(nn)[1]                       #number of individuals
  rm <- matrix(0,ncol=n,nrow=n)            #relatedness matrix
  for(i in 1:n){                        #probably a better way of doing this.
    tb <- table(nn[i,])
    rm[i,as.integer(names(tb))] <- tb/n
    ##     for(j in 1:n){
    ##       rm[i,j] <- mean(rm[i,]==j)
    ##     }
  }
  rm <- (rm+t(rm))/2                    #since nearest neighbour is not necessarily reflexive.
  diag(rm) <- 1
  return(rm)
}

#########################################################
##  kinship estimation 
#########################################################

kinship <- function(genotypes, method="correlation", ...){
  if(method=="emma"){
    K<-emma.kinship(t(genotypes), method="additive") 
  }
  else if(method=="correlation"){
    K <- cor(t(genotypes))
  }
  else if(method=="haploid"){
    K <- haploid.kinship(genotypes)
  }
  else if(method=="nearest.neighbour"){
    nn.matrix <- list(...)$nn.matrix
    K <- nn.relatedness(nn.matrix)
  }
  else{
    stop("Unknown method")
  }
  return(K)
}

#########################################################
## Mixed model analysis using the emma package
## http://mouse.cs.ucla.edu/emma/
#########################################################

mm.emma <- function(quant.trait, genotypes, kinship="correlation", ...){
  K<-kinship(genotypes, method=kinship, ... )
  mm<-emma.ML.LRT(t(quant.trait), t(genotypes), K)
  results <- matrix(mm$ps, ncol=1)
  return(results)
}

#########################################################
## Mixed model analysis using Matti Pirinen's mm
## package (needs to be available as ./mm/mm)
#########################################################

mm.mm <- function(quant.trait, genotypes, kinship="correlation", ...){

  K <- kinship
  if(!is.matrix(K)){
    K<-kinship(genotypes, method=kinship, ...)
  }

  n<-nrow(genotypes)
  nsnps <- ncol(genotypes)
  
  ## write out the files used for mm
  cat("GROUP1 ",n,"\n",sep="",file="data.mmtmp")
  write(t(quant.trait),ncol=1,file="data.mmtmp",append=TRUE)

  cat("GROUP1 ",n,"\n",sep="",file="cov.mmtmp")
  write(rep(1,n),ncol=1,file="cov.mmtmp",append=TRUE)

  cat("GROUP1"," ",n," ",1,"\n",sep="",file="K.mmtmp")
  write(t(K),ncol=n,file="K.mmtmp",append=TRUE)

  write.table(cbind(paste("SNP",seq(1,nsnps),sep=""),t(genotypes)),
              file="snps.mmtmp",col.names=FALSE,
              row.names=FALSE,quote=FALSE)

  ## run mm
  system('./mm/mm mm/params data.mmtmp cov.mmtmp snps.mmtmp K.mmtmp')
  system("sleep 2")                     #Just to allow file to sync.
  ## read results and cleanup
  read.table("snps.mmtmp_stats.out",as.is=TRUE,header=TRUE)->mm.res
  system("rm *.mmtmp")
  system("rm snps.mmtmp_stats.out")
  
  result <- mm.res[,c("chisq","snp_est")]
  result$chisq <- pchisq(result$chisq ,df=1,lower.tail=FALSE)
  names(result) <- NULL
  return(as.matrix(result))
}


#########################################################
## qq plots that compare corrected and uncorrected values
## 
#########################################################

qq.correction.plot <- function(uncorrected.results, corrected.results=list(), corrected.labels=c()){
  cols <- get.cols(length(corrected.results)+1)
  qq.plot.of.p.values(uncorrected.results[,1], col=cols[1], main="Corrections")

  c=2
  for(res in corrected.results){
    qq.plot.of.p.values(res[,1], col=cols[c], add=TRUE)
    c=c+1
  }
  legend("topleft", c("Uncorrected", corrected.labels), col=cols,pch=20,bty="n")
}

#########################################################
## collapse rare variants, either to a proportion, or to 
## an indicatior of one or more - typically use the former
## since it's more powerful according to Morris & Zeggini
## 2008.
## genotypes: nrow=individuals, ncol= variants. 
#########################################################

collapse.variants <- function(genotypes, method="proportion"){
  if(method=="proportion"){
    collapse <- apply(genotypes,1,mean)
  }
  else if(method=="indicator"){
    collapse <- apply(genotypes,1,max)
  }
  else{
    stop(paste("Unknown method", method))
  }
  return(collapse)
}


      



    

      
      
      
      
  
