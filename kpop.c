/*
  Copyright 2011 Iain Mathieson

  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  
  http://www.apache.org/licenses/LICENSE-2.0
  
  Unless required by applicable law or agreed to in writing, software 
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
  See the License for the specific language governing permissions and
  limitations under the License. 
*/

/* c implementation of the structured coalescent on a grid  */
/* This is basically a direct port of Gil McVean's R code */
/* Compile with: gcc -fpic -shared -o kpop.so kpop.c */
/* Load into R with dyn.load("kpop.so") and call using .C() */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)>(b)?(a):(b))

/* return a single sample from the integers 0:n-1 */
/* with probabilities given by probs need not sum to 1 */
int sample(int n, double *probs){
  int i;
  double cum_prob [n];
  double runif, p;

  p=0.0;
  for(i=0;i<n;++i){
    cum_prob[i]=p+probs[i];
    p=cum_prob[i];
  }

  runif=cum_prob[n-1]*rand()/(RAND_MAX+1.0);

  for(i=0;i<n;++i)
    if(runif<cum_prob[i])
      break;

  return i;
}

/* return the index of a randomly sample population */
/* which is equal to pop, but counting only the first*/
/* n_sum populations, uniformly sampled */
/* ==sample(which(pop.list==pop.for.migrant), 1) in R*/
int sample_populations(int *pop_list, int pop, int n_sum){
  int i,n, runif;
  int pop_indices [n_sum];
  
  n=0;
  for(i=0;i<n_sum;++i){
    if(pop==pop_list[i]){
      pop_indices[n]=i;
      n+=1;
    }
  }

  runif=floor(n*1.0*rand()/(RAND_MAX+1.0));
  return pop_indices[runif];
}

/* Structured coalescent simulation */
void simulate_tree_kpop_lattice_c(int *k_val, int *n, double *M, int *RNG_seed, double *Time, int *D1, int *D2, double* Branch_Length, int* Parent){

  int i,j,m, current_node, n_sum, pop_for_migrant, pop_new, pop_i, which_pop, pair_a, pair_b, cj;
  int k=*k_val;
  double t=0.0;
  int n_tot=0;
  int n_pop=k*k;
  int pop_ids [n_pop];
  double rate_co [n_pop];
  double rate_mig [n_pop];
  int this_mig_mat_row [4];
  double this_mig_prob_row [4];
  double rate_co_tot, rate_mig_tot, t_event, runif;

  cj=0;
  srand(*RNG_seed);

  for(i=0; i<n_pop; ++i){
    pop_ids[i]=i;
    n_tot+=n[i];
  }
  current_node=n_tot;
  n_sum=n_tot;

  int k_list [n_tot];
  for(i=0; i<n_tot; ++i){
    k_list[i]=i;
  }

  int pop_list [n_tot];
  m=0;
  for(i=0;i<n_pop;++i)
    for(j=0;j<n[i];++j)
      {
	pop_list[m]=i;
	++m;
      }

  if(*M<=0){
    printf("\n\n*** Error: Cannot have M<=0, setting to 0.01 ***\n\n");
    *M=0.01;
  };

  /* Store where lineages can migrate too */
  int mig_mat [n_pop][4];
  double mig_prob [n_pop][4];

  for(i=0;i<k;++i){
    for(j=0;j<k;++j){
      m=i*k+j;

      if(i>0){
	mig_mat[m][0]=(i-1)*k+j;
	mig_prob[m][0]=1;
      }
      else{
	mig_mat[m][0]=-1;
	mig_prob[m][0]=0.0;
      }

      if(i<k-1){
	mig_mat[m][1]=(i+1)*k+j;
	mig_prob[m][1]=1.0;
      }
      else{
	mig_mat[m][1]=-1;
	mig_prob[m][1]=0.0;
      }

      if(j>0){
	mig_mat[m][2]=i*k+j-1;
	mig_prob[m][2]=1.0;
      }
      else{
	mig_mat[m][2]=-1;
	mig_prob[m][2]=0.0;
      }

      if(j<k-1){
	mig_mat[m][3]=i*k+j+1;
	mig_prob[m][3]=1.0;
      }
      else{
	mig_mat[m][3]=-1;
	mig_prob[m][3]=0.0;
      }
    }
  }

  int mig_rates [n_pop];
  for(i=0;i<n_pop;++i){
    mig_rates[i]=(mig_mat[i][0]>=0)+(mig_mat[i][1]>=0)+(mig_mat[i][2]>=0)+(mig_mat[i][3]>=0);
    for(j=0;j<4;++j)
      mig_prob[i][j]=mig_prob[i][j]/mig_rates[i];
  }	

  while(n_sum>1){
    rate_co_tot=0.0;
    rate_mig_tot=0.0;
    for(i=0;i<n_pop;++i){
      rate_co[i]=n[i]*(n[i]-1.0)/2.0;
      rate_co_tot+=rate_co[i];
      rate_mig[i]=(*M/2.0) * n[i] * mig_rates[i];
      rate_mig_tot+=rate_mig[i];
    }
     
    runif=rand()/(RAND_MAX+1.0);
    t_event= -log(runif)/(rate_co_tot+rate_mig_tot);
    t+=t_event;

    /* Migration event */
    runif=rand()/(RAND_MAX+1.0);
    if(runif<=rate_mig_tot/(rate_co_tot+rate_mig_tot)){
      pop_for_migrant=sample(n_pop, rate_mig);
      pop_i=sample_populations(pop_list, pop_for_migrant, n_sum);

      for(j=0;j<4;++j){
	this_mig_mat_row[j]=mig_mat[pop_for_migrant][j];
	this_mig_prob_row[j]=mig_prob[pop_for_migrant][j];
      }
      pop_new=this_mig_mat_row[sample(4,this_mig_prob_row)];
      n[pop_list[pop_i]]-=1;
      n[pop_new]+=1;
      pop_list[pop_i]=pop_new;
    }
    /* Coalescent event */
    else{
      which_pop=sample(n_pop, rate_co);
      pair_a=sample_populations(pop_list, which_pop, n_sum);
      pair_b=sample_populations(pop_list, which_pop, n_sum);
      while(pair_b==pair_a){
	pair_b=sample_populations(pop_list, which_pop, n_sum);
      }

      Time[current_node]=t;
      /* R nodes are indexed from 1 */
      D1[current_node]=k_list[pair_a]+1;
      D2[current_node]=k_list[pair_b]+1;

      /* Branch lengths and Parent */
      Branch_Length[k_list[pair_a]]=t-Time[k_list[pair_a]];
      Branch_Length[k_list[pair_b]]=t-Time[k_list[pair_b]];
      Parent[k_list[pair_a]]=current_node+1;
      Parent[k_list[pair_b]]=current_node+1;

      k_list[MIN(pair_a,pair_b)]=current_node;
      k_list[MAX(pair_a,pair_b)]=k_list[n_sum-1];
      pop_list[MAX(pair_a,pair_b)]=pop_list[n_sum-1];
      
      n[which_pop]-=1;
      n_sum-=1;
      
      current_node+=1;
    }      
  }

  return;
}
  
/* Look at the daughters of node n. If they are terminal nodes, i.e. < n_sum, set the genotypes
to 1, otherwise keep going. Remember 1 based indexing! */
void fill_in_genotypes(int n_sum, int *gt, int *D1, int*D2, int node){
  int d1, d2;

  d1=D1[node]-1;
  d2=D2[node]-1;

  if(d1==-1 && d2==-1){		/* It's a child node */
    gt[node]=1;
  }
  else{
    if(d1<n_sum){
      gt[d1]=1;
    }
    else{
      fill_in_genotypes(n_sum, gt, D1, D2, d1);
    }
    
    if(d2<n_sum){
      gt[d2]=1;
    }
    else{
      fill_in_genotypes(n_sum, gt, D1, D2, d2);
    }
  }
}

/* Throw down a mutation on the tree and sets the genotypes to 1 for individuals with that mutation. 
Remember that for R compatability reasons the nodes are (yuk!) indexed from 0 */

void one_mutation_from_tree(int *n_sum, int *gt, int *D1, int *D2, double* branch_length){
  int this_mutation;
  int n_nodes=2*(*n_sum)-1;

  this_mutation=sample(n_nodes, branch_length);
  fill_in_genotypes(*n_sum, gt, D1, D2, this_mutation);
};

