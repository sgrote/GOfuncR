
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
// steffi:
#include "overall_sign.h" 

using namespace std ;

class go_groups {
	public:
		/**********
		 * save name and number of genes per group
		 ***********/
		go_groups( string &groups, string &sgenes, int co=1, string root_go = "GO:0003675" ) ;

		/**********
		 * calculates number of significant groups at different cutoffs
		 * returns array with these numbers of groups
		 * side effect: 1. saves the p-values to overall_significance
		 *              2. saves pvalues to data_pvals_{l,g}
		 ***********/
		int *calculate_data( string &data, double sum_nties, ostream *os=0 ) ;

		/**********
		 * equal to calculate_data, except side effect 2:
		 *              2. saves smallest pvalue over all groups to 
		 *					smallest_rand_p_{l,g}
		 ***********/
		int *calculate_rand( string &data, double sum_nties, ostream *os=0 ) ;

		/**********
		 * prints statistics to os, uses nr_randsets to calculate
		 * FWER. Runs overall_significance test and FDR estimate
		 * using less_sig and greater_sig.
		 ***********/
		void print_pvals( int nr_randsets, ostream &os ) ;
	private:
		vector<string> names ;
		vector<int> nr_of_genes ;
		// pvals of the dataset
		vector<double> data_pvals_l ;
		vector<double> data_pvals_g ;

		// smallest p-values of randomsets in both directions
		// needed for FWER, see print_pvals
 		multiset<double> smallest_rand_p_l ;
		multiset<double> smallest_rand_p_g ;

		int root_idx ;
		int cutoff ;
		
		overall_significance less_sig, greater_sig ;
	
} ;
