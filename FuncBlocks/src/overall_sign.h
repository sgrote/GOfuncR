
/* vim: set ts=4 tw=80 cindent : */

#ifndef OVERSIGN_H
#define OVERSIGN_H

#include <cstdlib>
#include <set>
#include <vector>
#include <map>
#include <iostream>

using namespace std ;


/****************************************************
  Methods to estimate significance
  from cumulative distribution functions of several
  datasets.
 */
class overall_significance {
	public: 
			/************
			 *  precision is the stepwidth for the pvalues.  
			 **************/
			overall_significance( double precision_=0.0001 ) ;

			/************
			 * saves a set to the private cdfs-vector.
			 * pvalues between [precision*m,precision*++m) will 
			 * be added up.
			 **************/
			void add_set( multiset<double> &pvals ) ;

			/************
			 * gives the overall significance of the cdfs-graph
			 * between 0 and co for dataset added as number set_.
			 **************/
			double significance( int set_, double co ) ;

			/************
			 * calculate the p-values at which the fdr is
			 * <0.1, <0.05, <0.01, <0.001 and <0.0001 and
			 * return these values. If it is not possible to 
			 * find a p-value between 0 and 1 that has a 
			 * fdr below the cutoff the value -1 will be returned.
			 **************/
			vector<double> *fdr( int set_ ) ;

			/*************
			 * calculate q values according to Yekutieli Benjamini 1999.
			 * see (9) p.179 
			 * returns vector of p-values. stepwidth according to precision.
			 ***************/
			// vector<double> *fdr_qvals( int set_ ) ;
			map<double,double> *fdr_qvals( int set_ ) ;

			/* for debugging */
			void print_cdfs( ostream &os ) ;
			double alt_sign( int ) ;

			// gives the index for the cdfs-vector for a pval 
			// index also works for fdr_qvals() returned vector
			int index_for_pval( double p ) ;
	private: 
			double precision ;
			int nr_sets ;
			
			// matrix of SETS x PVALUES 
			vector<vector<unsigned int>*> cdfs ; 
			vector<multiset<double> > cdfs_sets ;

} ;

#endif
