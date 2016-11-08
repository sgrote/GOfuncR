
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <set>
//#include "overall_sign.h"

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
		 * side effect: 1. saves the p-values to overall_significance
		 *              2. saves pvalues to data_pvals_{l,r}
		 ***********/
		int *calculate_data( ostream *os=0 ) ;

		/**********
		 * equal to calculate_data, except:
		 *   data is a string with random data from randomset
		 *   side effect 2:
		 *              2. saves smallest pvalue over all groups to 
		 *					smallest_rand_p_{l,r}
		 ***********/
		int *calculate_rand( string &data, ostream *os=0 ) ;

		/**********
		 * prints statistics to os, uses nr_randsets to calculate
		 * FWER. Runs overall_significance test and FDR estimate
		 * using osig_l and osig_r.
		 ***********/
		void print_pvals( int nr_randsets, ostream &os ) ;
	private:
		vector<string> names ; // GO ID
		vector<int> detected ; // # detected in group
		vector<int> changed_data ; // # changed in group

		// pvalues for all groups in dataset
		vector<double> data_pvals_l ;
		vector<double> data_pvals_r ;

		// smallest pvalue of each randomset
		multiset<double> smallest_rand_p_l ;
		multiset<double> smallest_rand_p_r ;

		//overall_significance osig_l, osig_r ;

		// index of root node
		int root_idx ;
		double cutoff ;
	
} ;
