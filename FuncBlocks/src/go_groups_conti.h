
#ifndef GO_GROUPS_CONTI_H
#define GO_GROUPS_CONTI_H

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include "overall_sign.h"

using namespace std ;

/*******************
 * encapsulates parsing and analysis of randomset-lines
 ********************/
class go_groups_conti {
	public:
		/**********
		 * save name and boolean (if more genes than co in group) for 
		 * each group
		 ***********/
		go_groups_conti( string &groups, istream *in=0, int co=1, string root = "GO:0003675" ) ;

		/**********
		 * calculates number of significant groups at different cutoffs
		 * returns array with these numbers of groups
		 * side effect: 1. saves the p-values to overall_significance
		 *              2. saves pvalues to data_pvals_{1,2}
		 ***********/
		int *calculate_data( string &data, ostream *os=0 ) ;

		/**********
		 * equal to calculate_data, except side effect 2:
		 *              2. saves smallest pvalue over all groups to 
		 *					smallest_rand_p_{1,2}
		 ***********/
		int *calculate_rand( string &data, ostream *os=0 ) ;

		/**********
		 * prints statistics to os, uses nr_randsets to calculate
		 * FWER. Runs overall_significance test and FDR estimate
		 * using osig_1 and osig_2.
		 ***********/
		void print_pvals( int nr_randsets, ostream &os ) ;
	private:
		vector<string> names ; // GO ID
		vector<bool> check ; // whether #genes in group > cutoff

		// pvalues for all groups in dataset
		vector<double> data_pvals_1 ;
		vector<double> data_pvals_2 ;

		// index of root node
		int root_idx ;

		// smallest pvalue of each randomset
		multiset<double> smallest_rand_p_1, smallest_rand_p_2 ;

		//overall_significance osig_1, osig_2 ;

} ;

#endif
