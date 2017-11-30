
#ifndef GO_GROUPS_HYPE_H
#define GO_GROUPS_HYPE_H

#include <string>
#include <vector>
#include <iostream>
#include <set>

using namespace std ;

/*******************
 * encapsulates parsing and analysis of randomset-lines
 ********************/
class go_groups_hyper {
	public:
		/**********
		 * save name, changed and detected 
		 ***********/
		go_groups_hyper( string &groups, string detected_s, string changed_s, string root_go = "GO:0003675", int cutoff_=1 ) ;

		/**********
		 * calculates number of significant groups at different cutoffs
		 * returns array with these numbers of groups
		 * side effect: 1. saves pvalues to data_pvals_{l,r}
		 ***********/
		int *calculate_data( ostream *os=0 ) ;

		/**********
		 * equal to calculate_data, except:
		 *   data is a string with random data from randomset
		 *   side effect: saves smallest pvalue over all groups to 
		 *					smallest_rand_p_{l,r}
		 ***********/
		int *calculate_rand( string &data, ostream *os=0 ) ;

		/**********
		 * prints statistics to os, uses nr_randsets to calculate
		 * FWER.
		 ***********/
		void print_pvals( int nr_randsets, ostream &os ) ;
	private:
		vector<string> names ; // GO ID
		vector<int> detected ; // # annotated to node
		vector<int> changed_data ; // # candidate annotated to node

		// NEW: expected and real number of annotated genes
		vector<double> n_anno_expected ;

		// pvalues for all groups in dataset
		vector<double> data_pvals_l ;
		vector<double> data_pvals_r ;

		// smallest pvalue of each randomset
		multiset<double> smallest_rand_p_l ;
		multiset<double> smallest_rand_p_r ;

		// index of root node
		int root_idx ;
		double cutoff ;
	
} ;

#endif
