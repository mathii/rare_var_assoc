#Functions to simulate genome-scale data with SMC and k by k populations on a lattice

local<-TRUE;
if (local) {
## 	setwd( "C:/Documents and Settings/mcvean/My Documents/Genealogies");
	source("coalescent.r");
} else {
	source("http://www.stats.ox.ac.uk/~mcvean/R/coalescent.r");
}

##########################################
#Function to find time in tree below node
##########################################

Tree.length<-function(tree, node, t=0) {

	if (tree[node,3]==0) {
		return(t);
	} else {
		t<-t+tree[node,2]-tree[tree[node,3],2];
		t<-tree.length(tree, tree[node,3], t);
		t<-t+tree[node,2]-tree[tree[node,4],2];
		t<-tree.length(tree, tree[node,4], t);
	}

	return(t);
}

##########################################
#2D color map for plotting trees
#returns an nxn map
##########################################

color.map.2d <- function(n){
  require(RColorBrewer)
##   x <- brewer.pal(n,"Spectral")
##   y <- rev(brewer.pal(n,"Spectral"))
##   x <- colorRampPalette(c("#E41A1C", "#377EBA"))(n)
  x <- colorRampPalette(c("#E41A1C", "#1A1CE4"))(n)
  y <- colorRampPalette(c("#000000", "#4DAF4A"))(n)
##   x <- colorRampPalette(c("#FF0000", "#0000FF"))(n)
##   y <- colorRampPalette(c("#000000", "#00FF00"))(n)

  x <- col2rgb(x)
  y <- col2rgb(y)

  fader <- abs(seq(-1,1,2/(n-1)))
  fader2 <- rev(seq(0,1,1/(n-1)))
  fader2 <- fader2*fader2
  fader3 <- sqrt(seq(0,1,1/(n-1)))
  
  n.cols<-matrix(0,nrow=3,ncol=n*n)
  for(i in 1:n){
    for(j in 1:n){
      n.cols[,n*(i-1)+j] <- pmax((1-fader2[j])*fader[i]*x[,i]+(fader2[j])*x[,i],fader3[j]*y[,j])
    }
  }

  return(mapply(rgb, n.cols[1,], n.cols[2,], n.cols[3,], MoreArgs=list(maxColorValue=255)))
}

#################################################################
#Function to plot a tree that comes with migration information
#Note tree HAS to be time ordered
#################################################################

plot.tree<-function(data, maxdepth=-1, add.mutations=FALSE, y.low=0, names.plot=TRUE, lwd=1, cex.mtn=0.5, 
	plot.order=c(),cex.txt=0.5, srt.txt=90, annot=c(), n.pop=1, col="black", return.plot.order=FALSE,
	add.internal=TRUE) {

	nseq<-(nrow(data)+1)/2;
	mrca = nrow(data);
	if (length(plot.order)==0) plot.order<-order.seqs(data, node=mrca);
	is.mig=FALSE;

	#Choose colours for popns
	if (n.pop>1 & length(annot)>0) {
		is.mig<-TRUE;
## 		pop.pal<-sample(rainbow(n.pop));
                pop.pal <- color.map.2d(sqrt(n.pop))
              }

	if (maxdepth<0) maxdepth=max(data[,2]);
	plot(0,0, type="n", xaxt="n", xlab="", ylab="Time", bty="n", ylim=c(y.low,maxdepth), xlim=c(0,nseq+1));

	#Get positions of terminal branches
	branch.pos<-order(plot.order);

	#Iterate through nodes finding non-terminal ones
	for (node in (nseq+1):mrca) if (data[node,3]>0) {

		#Find branch pos
		branch.pos[node]<-(branch.pos[data[node,3]]+branch.pos[data[node,4]])/2;

		if (add.internal) text(node, x=branch.pos[node], y=data[node,2], cex=cex.txt, pos=1);

		#DO horizontal line
		if (is.mig) col<-pop.pal[as.integer(annot[[node]][1,1])];
		segments(x0=branch.pos[data[node,3]], y0=data[node,2], x1=branch.pos[data[node,4]], y1=data[node, 2], 
			lwd=lwd, col=col);

		#DO LH vertical
		x0<-branch.pos[data[node,3]];
		d1<-data[node,3];
		segments(x0, data[d1,2], x0, data[node, 2], col=col, lwd=lwd);

	#If migrations, need to add in descendants lines
		if (is.mig) for (j in 1:ncol(annot[[d1]])) {
			col<-pop.pal[as.integer(annot[[d1]][1,j])];
			segments(x0, annot[[d1]][2,j], x0, data[node,2], col=col, lwd=lwd);
		}

		#Add on mutations if wanted
		if (add.mutations==T && data[data[node,3],5]>0) {
			time.span<-c(data[data[node,3],2], data[node,2]);
			y.vals<-time.span[1]+runif(data[data[node,3],5])*(time.span[2]-time.span[1]);
			points(rep(branch.pos[data[node,3]], data[data[node,3],5]), y.vals, pch=19, col=rgb(0,0,0.5), cex=cex.mtn);
		}

		#DO RH vertical
		x0<-branch.pos[data[node,4]];
		d2<-data[node,4];
		segments(x0, data[d2,2], x0, data[node, 2], col=col, lwd=lwd);
	
		#If migrations, need to add in descendants lines
		if (is.mig) for (j in 1:ncol(annot[[d2]])) {
			col<-pop.pal[as.integer(annot[[d2]][1,j])]
			segments(x0, annot[[d2]][2,j], x0, data[node,2], col=col, lwd=lwd);
		}

		#Add in mutations if wanted
		if (add.mutations==T && data[data[node,4],5]>0) {
			time.span<-c(data[data[node,4],2], data[node,2]);
			y.vals<-time.span[1]+runif(data[data[node,4],5])*(time.span[2]-time.span[1]);
			points(rep(branch.pos[data[node,4]], data[data[node,4],5]), y.vals, pch=19, col=rgb(0,0,0.5), cex=cex.mtn);
		}

		#End of loop for node
	}

	#Add names of sequences at bottom
	if (names.plot==TRUE) {
		text(x=c(1:nseq), y=rep(0, nseq), labels=data[plot.order,1], srt=srt.txt, cex=cex.txt, pos=1)
	}

	if (return.plot.order) return(plot.order)
}

################
#Arrangement
#1 2 3 ... k
#(k+1) (k+2)...
#.........k^2
################

####################################
#Each pop has and i,j location
#Also an identifier m = (i-1)*k+j;
#This function decodes m;
####################################

decode.m<-function(m, k=2) {
	if (m>k^2) {
		cat("\n\n***Error: index is larger than max possible value***\n\n");
		return(NULL);
	}
	if (m<1) {
		cat("\n\n***Error: index is <1***\n\n");
		return(NULL);
	}
	i<-as.integer((m-1)/k)+1;
	j<-m-(i-1)*k;
	return(c(i,j));
}

#########################################################
# For each individual, compute their nearest neighbour
# (i.e. first coalescence) from a tree object . Choose
# randomly when there are ties. 
#########################################################

nearest.neighbour.vector <- function(t){
  n <- (dim(t$tree)[1]+1)/2
  t$tree <- cbind(t$tree, closest=-1, distance=-1)

  for(i in (n+1):(2*n-1)){
    D1=t$tree[i,"D1"]
    D2=t$tree[i,"D2"]
    if(D1<=n & D2<=n){
      t$tree[D1,"closest"] <- D2
      t$tree[D1,"distance"] <- 0
      t$tree[D2,"closest"] <- D1
      t$tree[D2,"distance"] <- 0
      t$tree[i,"closest"] <- sample(c(D1,D2),1)
      t$tree[i,"distance"] <- 1
    }
    else if(D1<=n & D2>n){
      t$tree[D1,"closest"] <- t$tree[D2,"closest"]
      t$tree[D1,"distance"] <- 0
      t$tree[i,"closest"] <- D1
      t$tree[i,"distance"] <- 1
    }
    else if(D2<=n & D1>n){
      t$tree[D2,"closest"] <- t$tree[D1,"closest"]
      t$tree[D2,"distance"] <- 0
      t$tree[i,"closest"] <- D2
      t$tree[i,"distance"] <- 1
    }
    else if(D1>n & D2>n){
      new.closest <- -1
      new.distance <- -1
      if(t$tree[D1,"distance"]>t$tree[D2,"distance"]){
        idx <- D1
      }
      else if(t$tree[D1,"distance"]<t$tree[D2,"distance"]){
        idx <- D2
      }
      else if(t$tree[D1,"distance"]==t$tree[D2,"distance"]){
        idx=sample(c(D1,D2),1)
      }
      else{
        stop("Shouldn't be here")
      }

      t$tree[i,"closest"] <- t$tree[idx,"closest"]
      t$tree[i,"distance"] <- t$tree[idx,"distance"]+1
    } 
    else{
      stop("Shouldn't be here")
    } 
  }

  return(t$tree[1:n,"closest"])
}
      

##################################################
## Simulate samples on a k by k population lattice 
##################################################

simulate.tree.kpop.lattice<-function(k=2, n=rep(10,k^2), M=1) {

	#cat("\n\nSimulating a tree for a population lattice\n\n");

	if (M==0) {
		cat("\n\n*** Error: Cannot have M=0 ***\n\n", sep="");
	}
	if (M>10) {
		cat("\n\n*** Warning: M>10, resetting M=10 ***\n\n", sep="");
                M=10;
	}

	t<-0;
	n.pop<-k^2;
	pop.ids<-1:n.pop;
	n.tot<-sum(n);
	tree<-matrix(0, nrow=2*n.tot-1, ncol=5);
	colnames(tree)<-c("Node.ID", "Time", "D1", "D2", "Mutations");
	tree[,1]<-c(1:nrow(tree));
	current.node<-n.tot+1;
	k.list<-c(1:n.tot);
	pop.list<-rep(1:n.pop, times=n);  

	#Each population m (1 to n.pop) has an i,j entry - see above function        
	#Store where lineages can migrate to
	mig.mat<-array(0, c(n.pop,4));
	for (i in 1:k) for (j in 1:k) {
		m<-(i-1)*k+j;
		if (i>1) mig.mat[m,1]<-(i-2)*k+j;
		if (i<k) mig.mat[m,2]<-i*k+j;
		if (j>1) mig.mat[m,3]<-(i-1)*k+j-1;
		if (j<k) mig.mat[m,4]<-(i-1)*k+j+1;
	}

	mig.rates<-apply(mig.mat>0, 1, sum);
	mig.prob<-(mig.mat>0);
	mig.prob<-mig.prob/apply(mig.prob, 1, sum);

	#Initialise information about location
	annot<-list();
	for (i in 1:(2*n.tot-1)) {
		annot[[i]]<-array(0, c(2,1));
		if (i<=n.tot) annot[[i]][,1]<-c(pop.list[i], 0);
	}

        cj=0
        ct=0

	while(sum(n)>1) {

		#Work out rate of events and choose time for next
		rate.co<-n*(n-1)/2;
                
		rate.co.tot<-sum(rate.co);
		rate.mig<-M/2*n*mig.rates;
		rate.mig.tot<-sum(rate.mig);
		t.event<--log(runif(1))/(rate.co.tot+rate.mig.tot);
                
		t<-t+t.event;

##                 if(sum(n)==2){
## ##                   print(pop.list)
##                               cj=cj+1
##                               ct=ct+t.event
##                               print( c(pop.list, rate.co.tot) )
##                             }

                
		#Migration event
		if ((runif(1)<=rate.mig.tot/(rate.mig.tot+rate.co.tot))){
                  #Find popn where migrant occurs
                        pop.for.migrant<-sample(pop.ids, 1, prob=rate.mig);
                        sample.from <- which(pop.list==pop.for.migrant)
                        if(length(sample.from)>1){ i<-sample(sample.from, 1) }
                        else{i<-sample.from}
        	        pop.new<-sample(mig.mat[pop.for.migrant,], 1, prob=mig.prob[pop.for.migrant,]);
			n[pop.list[i]]<-n[pop.list[i]]-1;
			n[pop.new]<-n[pop.new]+1;
			pop.list[i]<-pop.new;
			annot[[k.list[i]]]<-cbind(annot[[k.list[i]]], c(pop.new, t));
		}
		
		#Coalescent event
		else {
			#cat("\rDown to ", sum(n)-1, " lineages", sep="");
			which.pop<-sample(pop.ids, 1, prob=rate.co); 
			pair<-sample(which(pop.list==which.pop), 2);
			tree[current.node,2]<-t;
			tree[current.node,3]<-k.list[pair[1]];
			tree[current.node,4]<-k.list[pair[2]];
			k.list[min(pair)]<-current.node;
			k.list[max(pair)]<-k.list[length(k.list)];

			pop.list[max(pair)]<-pop.list[length(pop.list)];
			n[which.pop]<-n[which.pop]-1;

			annot[[current.node]][,1]<-c(which.pop, t);
	
			current.node<-current.node+1;

			k.list<-k.list[1:sum(n)];
			pop.list<-pop.list[1:sum(n)];

		}
	}

	#cat("\n\nDone!\n\n");
##         print(c(cj, ct))
	return(list(tree=tree, annot=annot));	

}

#########################################
#Same as above but tries to call the c
#implementation of the above ( in kpop.c )
#This actually can do the annotation in
#annotate.tree, but doesn't have the annot
#in the list like the non c function
#########################################

simulate.tree.kpop.lattice.c<-function(k=2, n=rep(10,k^2), M=1, annotate.tree=FALSE) {
  nr=2*sum(n)-1
  Time=as.double(rep(0,nr))
  D1=as.integer(rep(0,nr))
  D2=as.integer(rep(0,nr))
  Branch.length=as.double(rep(0,nr))
  Parent=as.integer(rep(0,nr))
  
  seed=as.integer(floor(runif(1)*65535))
  
  res<-.C("simulate_tree_kpop_lattice_c",as.integer(k), as.integer(n), as.double(M), seed, Time, D1, D2, Branch.length, Parent )

  if(annotate.tree){
    result<-cbind(Node.ID=1:nr, Time=res[[5]], D1=res[[6]], D2=res[[7]], Mutations=0, Branch.length=res[[8]], Parent=res[[9]])
  }
  else{
    result<-cbind(Node.ID=1:nr, Time=res[[5]], D1=res[[6]], D2=res[[7]], Mutations=0)
  }
  
  return(list(tree=result))
}
  
#########################################
#Function to add branch lengths to tree
#########################################

annotate.tree<-function(tree) {
	colnames.old<-colnames(tree);
	tree<-cbind(tree, array(0, c(nrow(tree), 2)));
	colnames(tree)<-c(colnames.old, "Branch.length", "Parent");
	which.add<-ncol(tree)-1;
	which.add.2<-which.add+1;
	n<-(nrow(tree)+1)/2;
	for (j in (n+1):(2*n-1)) {
		tree[tree[j,3],which.add]<-tree[j,2]-tree[tree[j,3],2];
		tree[tree[j,4],which.add]<-tree[j,2]-tree[tree[j,4],2];
		tree[tree[j,3],which.add.2]<-j;
		tree[tree[j,4],which.add.2]<-j;
	}
	return(tree);
}



###########################################################
#Function to list lineages present at a given time point
###########################################################

lineages.time<-function(tree, time) {

	mrca<-nrow(tree);

	op<-c();
	for (i in mrca:(nseq+1)) {
		if (tree[i,2]>time) {
			if (tree[tree[i,3],2]<time) op<-c(op, tree[i,3]);
			if (tree[tree[i,4],2]<time) op<-c(op, tree[i,4]);
		}
	}
	return(op);
}

######################################################################
#Function to take a tree and choose a point for recombination event
#written in simple way so that generalised version could be included
######################################################################

sample.rec.event<-function(tree, annot=c(), branch.length.col=6) {

	#Choose branch for recombination event
	which.branch<-sample(1:nrow(tree), 1, prob=tree[,branch.length.col]);

	return(which.branch);
}

##################################################################################
#Function to sample time for recombination event on a branch in an annotated tree
#Could generalise to allow inhomogenous rec rates (e.g. due to varying Ne);
##################################################################################

choose.rec.time<-function(tree, annot=c(), which.branch, branch.length.col=6, parent.col=7) {
	t.rec<-tree[which.branch,2]+runif(1)*(tree[tree[which.branch,parent.col],2]-tree[which.branch,2]);
	return(t.rec);
}

##################################################################################
#Function to sample a coalescent event for a floating lineage
##################################################################################

choose.co.event<-function(tree, which.branch, t.rec,  annot=c(), branch.length.col=6, parent.col=7, debug=TRUE) {

	#Number of sequences

	#Get list of lineages present at recombination event
	lins<-lineages.time(tree, t.rec);

	#Remove recombined lineage
	lins<-lins[lins!=which.branch];

	n.lin<-length(lins);  #-1 because must include recombined lineage
	epoch<-2*nseq-length(lins);
	is.free<-TRUE;
	t.curr<-t.rec;

	if (debug) cat("\n\n*** Sampling coalescent event for lineage ", which.branch, " recombining at ", t.rec, " ***\n", sep="");

	#Now iterate up tree - this version only considers 1 popn
	while(is.free) {

		#Choose time for next event
		t.co<- -log(runif(1))/n.lin;
		t.curr<-t.co+t.curr;
	
		#If a co event - choose who to coalesce with (be careful with 'sample' function).
		if (t.curr<tree[epoch,2]) {
			is.free<-FALSE;
			if (n.lin>1) {
				who.co<-sample(lins, 1);
			} else {
				who.co<-lins[1];
			}
			if (debug) cat("\nFree lineage ", which.branch, " coalesces with ", who.co, " at ", t.curr, sep="");

		#Else if time>next change in # of lineages - restart
		} else {
			t.curr<-tree[epoch,2];

			#remove daughter lineages from lins
			lins<-setdiff(lins, tree[epoch,3:4]);
			lins<-c(lins, epoch);
			n.lin<-length(lins);

			if (debug) cat("\nNext event is a coalescent in tree: ", tree[epoch,3], "^", tree[epoch,4], " at ", t.curr, sep="");

			epoch<-epoch+1;
		}

		#If has gone beyond MRCA just sample;
		if (epoch>mrca) {
			t.co<- -log(runif(1));
			t.curr<-t.curr+t.co;
			is.free<-FALSE;
			who.co<-mrca;

			if (debug) cat("\n**Reached MRCA in tree - sampling time to new MRCA", sep="");
		}
	}

	if (debug) cat("\n..Finished\n");
	return(list(who.co=who.co, t.co=t.curr));
}

#####################################################
#Function to prune an re-graft a branch from a tree
#####################################################

prune.regraft<-function(tree, which.branch, t.rec, who.co, t.co, debug=TRUE) {

	nseq<-(nrow(tree)+1)/2;
	mrca<-nrow(tree);

	if (debug) cat("\n\n** Prune-regraft operation to move node ", which.branch, " at ", 
		t.rec, " to co with ", who.co, " at ", t.co, "\n\n", sep="");

	#The following could be used to prevent unnecessary operations under certain conditions
	#If rec in top of tree, very simple operation - change MRCA;	
	#If coalesce with parental node, also simple
	#If coalesce with sister node, also simple

	#Find parent, grandparent and sister nodes - note if pn==mrca, gpn=pn
	pn<-tree[which.branch,7];
	dd1<-as.integer(tree[pn,4]==which.branch);
	sister.node<-tree[pn,4-dd1];
	if (pn!=mrca) {
		gpn<-tree[pn,7];
	} else {
		gpn<-pn;
	}
	dd2<-as.integer(tree[gpn,4]==pn);

	if (debug) cat("PN = ", pn, ": GPN = ", gpn, ": sister = ", sister.node, "\n", sep="");

	#Now find parental and sister nodes for branch coalesced with
	pn.co<-tree[who.co,7];

	#Note if coalesce with sister node, their parent is actually your grandparent
	if (who.co==sister.node) pn.co<-gpn;
	if (who.co!=mrca) {
		dd3<-as.integer(tree[pn.co,4]==who.co);
		if (who.co==sister.node) dd3<-as.integer(tree[pn.co,4]==pn);
	}

	if (debug) cat("Parent of co = ", pn.co, "\n", sep="");

	#Currently make a copy, but no need ultimately	
	tree2<-tree;

	#Now revise nodes in tree - new node will replace pn
	#First remove parental node
	tree2[gpn,3+dd2]<-sister.node;
	tree2[sister.node,7]<-gpn;
	tree2[sister.node,6]<-tree2[gpn,2]-tree2[sister.node,2];

	#If you are coalescing with your parental node, you are coalescing with your sister node;
	if (pn==who.co) who.co<-sister.node;

	if (debug) cat("GPN node ", gpn, " gets sister ", sister.node, " as new daughter\n");

	#Next insert new node below parent of coalescent partner
	tree2[who.co,7]<-pn;
	if (who.co!=mrca) tree2[pn.co,3+dd3]<-pn;

	if (debug) cat("PN node of co ", pn.co, " gets new daughter ", pn, "\n", sep="");

	#Now make new node
	tree2[pn,2]<-t.co;
	tree2[pn,3]<-which.branch;
	tree2[pn,4]<-who.co;
	tree2[pn,7]<-pn.co;
	if (pn.co>0) {
		tree2[pn,6]<-tree2[pn.co,2]-tree2[pn,2];
	} else {
		tree2[pn,6]<-0.0;
	}
	
	tree2[which.branch,7]<-pn;
	tree2[who.co,7]<-pn;
	tree2[which.branch,6]<-tree2[pn,2]-tree2[which.branch,2];
	tree2[who.co,6]<-tree2[pn,2]-tree2[who.co,2];

	#Now need to sort nodes in tree
	#This is in efficient - there is actually a simpler structure that can be exploited
	#Basically need to subtract 1 off nodes that lie between t.rec and t.co

	tree2<-tree2[order(tree2[,2]),];
	mm<-match(1:mrca, tree2[,1]);
	ii<-(nseq+1):mrca;
	tree2[ii,3]<-mm[tree2[ii,3]];
	tree2[ii,4]<-mm[tree2[ii,4]];
	tree2[1:(mrca-1),7]<-mm[tree2[1:(mrca-1),7]];

	tree2[mrca,6]<-0;
	tree2[mrca,7]<-0;

	tree2[,1]<-1:mrca;
	
	return(tree2);
}




#Function to take a tree and simulate recombination event
simulate.next.tree.smc<-function(tree, annot=c(), branch.length.col=6, show.tree=TRUE, parent.col=7) {

	#Number of sequences
	nseq<-(nrow(tree)+1)/2;
	mrca<-nrow(tree);

	#Choose branch for recombination event
	which.branch<-sample.rec.event(tree, annot, branch.length.col);

	#Choose time for recombination event
	t.rec<-choose.rec.time(tree, annot, which.branch, branch.length.col, parent.col);

	#Choose time and node for coalescent event
	co.event<-choose.co.event(tree, which.branch, t.rec, annot, branch.length.col, parent.col);

	#Prune and re-graft tree
	tree2<-prune.regraft(tree, which.branch, t.rec, co.event$who.co, co.event$t.co);

	if (show.tree) plot.tree(tree2);
	
	return(tree2);
}



#Function to check tree...

check.tree<-function(tree) {

	#Number of sequences
	nseq<-(nrow(tree)+1)/2;
	mrca<-nrow(tree);
	fl<-0;

	#Check that times are monotonic
	if (min(diff(tree[,2]))<0) {
		cat("\nError: times are not monotonic\n");
		f1<-1;
	}

	#Check ancestors and descendants;
	for (i in (nseq+1):mrca) {
		if (tree[tree[i,3],7]!=i | tree[tree[i,4],7]!=i) {
			cat("\nError: descendants do not tally\n");
			fl<-1;
		}
	}
	
	#Check branch lengths
	for (i in 1:(mrca-1)) {
		tl<-tree[tree[i,7],2]-tree[i,2];
		if (abs(tl-tree[i,6])>1e-6) {
			cat("\nError: branch lengths do not agree for node ", i, "\n");
			fl<-1;
		}
	}

	return(fl);
}

############################################################
#Function to choose a location on a tree to add a mutation
#The uniform argument is confusing. If true, it means that
#mutations occur uniformly along each branch, so you get
#a 1/f frequency spectrum. If false, mutations occur at
# each node proportional to its age, so you get a more
# uniform distribution.
############################################################

choose.mutation<-function(tree, uniform=TRUE, n.mutations=1) {

  nr <- nrow(tree)
  if (uniform){return(sample(1:nr, n.mutations, prob=tree[,6], replace=TRUE))}
  else{
    pr <- sapply(t$annot,dim)[2,]
    pr <- 2^pr
    pr[nr]<-0
    return(sample(1:nr, n.mutations, prob=pr, replace=TRUE))
  }
}  

##############################################################################
#Function to add a mutation to a given node and return the daughter sequences
##############################################################################

add.mutation<-function(tree, node) {
	n<-(nrow(tree)+1)/2;
	seq<-rep(0, n);
	seq[find.descendants(node, tree, c())]<-1;
	return(seq);
}

####################################################
#Funciton to plot spatial distribution of mutations
####################################################

plot.mutation<-function(seq, k=1, n=c(length(seq)), n.bks=10, seg.col="black", override.pal=NULL, ...) {
  plot(0,0,type="n",xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
  stamp.mutation(seq, k, n, n.bks, seg.col, override.pal=override.pal)
}

####################################################
#Funciton to plot spatial distribution of mutations
#on top of another plot. 
####################################################

stamp.mutation<-function(seq, k=1, n=c(length(seq)), n.bks=10, seg.col="black", x.scale=1, x.shift=0, y.scale=1, y.shift=0, override.pal=NULL, ...) {

  if(require(RColorBrewer)){
    n.bks <- min(n.bks, 9)
    pal<-brewer.pal(n.bks,"YlOrRd")
  }
  else{
    pal<-heat.colors(n.bks);
                                        #pal<-grey((1:n.bks)/(n.bks+1)));
    pal<-pal[length(pal):1];
  }
  pal <- c(pal, "black")
  
	pos<-(0:k)/k;
	bks<-seq(-0.01, 1.01, length=n.bks+1);

	rect(0+x.shift,0+y.shift,1*x.scale+x.shift,1*y.scale+y.shift, col=seg.col);
	if (k>1) {
		for (i in 1:(k-1)) {
			segments((i/k)*x.scale+x.shift, 0+y.shift, (i/k)*x.scale+x.shift, 1*y.scale+y.shift, col=seg.col);
			segments(0+x.shift,(i/k)*y.scale+y.shift,1*x.scale+x.shift,(i/k)*y.scale+y.shift, col=seg.col);
		}
	}

	l.freq<-aggregate(seq, list(rep(1:(k^2), each=n)), mean);
	for (i in 1:k^2) {
		id<-decode.m(i, k);
		if(any(is.null(override.pal))){
                  rect(pos[id[1]]*x.scale+x.shift,
                       pos[id[2]]*y.scale+y.shift,
                       pos[id[1]+1]*x.scale+x.shift,
                       pos[id[2]+1]*y.scale+y.shift,
                       col=pal[findInterval(l.freq[i,2], bks)], border=seg.col);
                }
                else{
                  rect(pos[id[1]]*x.scale+x.shift,
                       pos[id[2]]*y.scale+y.shift,
                       pos[id[1]+1]*x.scale+x.shift,
                       pos[id[2]+1]*y.scale+y.shift,
                       col=override.pal[i], border=seg.col);
                } 
              }
}


#################################################
#Function to summarise sequences over a lattice
#################################################

summarise.lattice<-function(seq, k=1, n=length(seq)) {

	#Fst
	l.freq<-aggregate(seq, list(rep(1:(k^2), each=n)), mean);
	mn<-mean(seq);
	n.t<-length(seq);
	pwd.w<-mean(2*l.freq[,2]*(1-l.freq[,2])*n/(n-1));
	pwd.t<-2*mn*(1-mn)*n.t/(n.t-1);
	fst<-1-pwd.w/pwd.t;

	return(list(x.bar=sum(seq), fst=fst));

}

#################################################
#Sample individuals, get n.tot individuals
#from a k by k grid and return a vector
#of length k*k with the number of individuals
#in each box. 
#################################################

sample.individuals <- function(k, n.tot){
  positions<-sample(k*k, n.tot, replace=TRUE)
  tab<-table(positions)
  counts<-rep(0,k*k)
  counts[as.numeric(names(tab))]<-tab
  return(counts)
}

#################################################
#Decode counts. For a list of k*k cell counts
#return the position of each element
#################################################

decode.cell.counts <- function(k, n.sim){
  counts <- sapply( 1:(k*k), decode.m, k=k ) 
  x<-rep(counts[1,], times=n.sim)
  y<-rep(counts[2,], times=n.sim)
  return(cbind(x,y))
}
  
## set.seed(234343);

## k.sim<-5;
## n.sim<-rep(5, k.sim^2);
## #n.sim<-c(5,0,0,0);
## t<-simulate.tree.kpop.lattice(k=k.sim, n=n.sim, M=0.5);

## dev.new()
## plot.tree(t$tree, annot=t$annot, n.pop=k.sim^2, lwd=2, add.internal=F);

## #Annotate tree with branch lengths and parent nodes
## t$tree<-annotate.tree(t$tree);

## #Plot
## dev.new()
## plot.mutation(add.mutation(t$tree, choose.mutation(t$tree)), k.sim, n.sim);


## #Series of simulations
## k.sim<-5;
## n.sim<-rep(5, k.sim^2);
## n.s<-100;
## k.s<-10;
## op<-array(0, c(n.s*k.s, 2));
## colnames(op)<-c("Freq", "Fst");
## ct<-0;
## cat("\n\nSimulations\n\n");
## for (i in 1:n.s) {
## 	cat("\nTree ", i, "\t", sep="");
## 	t<-simulate.tree.kpop.lattice(k=k.sim, n=n.sim, M=10);
## 	t$tree<-annotate.tree(t$tree);
## 	for (j in 1:k.s) {
## 		ct<-ct+1;
## 		seq<-add.mutation(t$tree, choose.mutation(t$tree));
## 		tmp<-summarise.lattice(seq, k.sim, n.sim);
## 		op[ct,1]<-tmp$x.bar;
## 		op[ct,2]<-tmp$fst;		
## 	}
## }
## cat("\n\nDone!\n\n");

## dev.new()
## par(mfrow=c(1,2));

## neu<-1/(1:(sum(n.sim)-1));
## neu<-neu/sum(neu);

## nn<-sample(1:(sum(n.sim)-1), nrow(op), replace=T, prob=neu);
## qqplot(nn, op[,1], xlab="Neutral", ylab="Lattice");
## abline(0,1,lty="dotted", col="red");

## plot(op);
## bks<-20;
## bks<-seq(-0.01, 1.01, length=bks+1);
## mn.fst<-aggregate(op[,2], list(findInterval(op[,1]/sum(n.sim), bks)), mean);
## mn.f<-aggregate(op[,1], list(findInterval(op[,1]/sum(n.sim), bks)), mean);
## points(mn.f[,2], mn.fst[,2], pch=19, col="red", type="b");

## par(mfrow=c(1,1));

## p<-0.2;
## op.chk<-c();
## for (i in 1:10) {
## 	tmp<-rbinom(sum(n.sim), 1, p);
## 	xx<-summarise.lattice(tmp, k.sim, n.sim);
## 	op.chk[i]<-xx$fst;
## }

## n.s<-100;
## tree<-t$tree;
## for (i in 1:n.s) {
## 	cat("\rIteration ", i, "     ", sep="");
## 	tree2<-simulate.next.tree.smc(tree);
## 	fl<-check.tree(tree2);
## 	if (fl) break();
## 	Sys.sleep(0.1);
## 	tree<-tree2;
## }




