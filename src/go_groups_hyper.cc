
#include "go_groups_hyper.h"
#include <cmath>
#include <map> 
#include <set> 
#include <fstream>
#include <iomanip>      // std::setprecision

#include <Rcpp.h>

using namespace std ;

go_groups_hyper::go_groups_hyper( string &groups, string detected_s, string changed_s, string root_go, int cutoff_ )  
{
	cutoff = static_cast<double>(cutoff_) ;
	istringstream is( groups.c_str() ) ;
	istringstream is_det( detected_s.c_str() ) ;
	istringstream is_ch( changed_s.c_str() ) ;
	string name ;
	int det,ch ;
	//Rcpp::Rcout << endl << "go_groups_hyper:" << endl;
	while ( is >> name ) {
		is_det >> det ;
		is_ch >> ch ;		
		//Rcpp::Rcout << name << "\t" << det << "\t" << ch << endl ;
		names.push_back( name ) ;	// vector<string> names: node-IDs
		detected.push_back( det ) ; // vector<int> detected: # genes annotated to node
		changed_data.push_back( ch ) ; // vector<int> changed_data: # candidate genes in node
		if ( name == root_go ) root_idx = names.size() - 1 ; // index of root node
	}
	//Rcpp::Rcout << "GOs: " << names.size() << endl ;
//	Rcpp::Rcerr << detected.size() << endl ;
//	Rcpp::Rcerr << changed_data.size() << endl ;


}


int* go_groups_hyper::calculate_data( ostream *os ) 
{
	// count number of nodes with p-value < 0.1, 0.05, ... 
	int *ret = new int[10] ;
	for ( int i=0 ; i<10 ; ++i ) {
		ret[i] = 0 ;
	}

	vector<int> &changed = changed_data ;

	data_pvals_l.resize( names.size() ) ; // vector<double> data_pvals_l: p-vals for all nodes 
	data_pvals_r.resize( names.size() ) ;
	
	n_anno_expected.resize( names.size() ) ; // vector<double>: expected no. of anno. candidate genes per GOs 

	//multiset<double> pvals_l, pvals_r ; // set of p-vals

	// loop over all nodes, get p-value for over- and underrep
	for( unsigned int i = 0 ; i < names.size() ; ++i ) {
		
		data_pvals_l[i] = 2. ;
		data_pvals_r[i] = 2. ;

//		Rcpp::Rcerr << i << endl ;

//		Rcpp::Rcerr << detected.size() << endl ;

		double N = detected[root_idx] ;
		double M = changed[root_idx] ;
		double n = detected[i] ;
		double x = changed[i] ;
		

		if ( n < cutoff ) continue ; 

		double prob_left;
		double prob_right ;

		// pval
		prob_left  = R::phyper( x, M, N-M, n, 1, 0 ) ; // underrep
		prob_right = R::phyper( x-1., M, N-M, n, 0, 0 ) ; // overrep
		//prob_left  = phyper( x, M, N-M, n, 1, 0 ) ; 
		//prob_right = phyper( x-1., M, N-M, n, 0, 0 ) ;
		data_pvals_l[i] = prob_left ;
		data_pvals_r[i] = prob_right ;

		//pvals_l.insert( prob_left ) ;
		//pvals_r.insert( prob_right ) ;
		
		// NEW: expected no. of candidate genes
		n_anno_expected[i] = n*(M/N);

		//if ( os ) {
			//*os << names[i] << "\t" 
			   //<< N << "\t"
			   //<< n << "\t"
			   //<< M << "\t"
			   //<< x << "\t" << endl ;
		//}
		// count number of nodes with p-value < 0.1, 0.05, ... 
		if ( prob_left < 0.1 ) {
			ret[0]++ ;
			if ( prob_left < 0.05 ) {
				ret[1]++ ;
				if ( prob_left < 0.01 ) {
					ret[2]++ ;
					if ( prob_left < 0.001 ) {
						ret[3]++ ;
						if ( prob_left < 0.0001 ) {
							ret[4]++ ;
						}
					}
				}
			}
		}
		if ( prob_right < 0.1 ) {
			ret[5]++ ;
			if ( prob_right < 0.05 ) {
				ret[6]++ ;
				if ( prob_right < 0.01 ) {
					ret[7]++ ;
					if ( prob_right < 0.001 ) {
						ret[8]++ ;
						if ( prob_right < 0.0001 ) {
							ret[9]++ ;
						}
					}
				}
			}
		}
	}
	//osig_l.add_set( pvals_l ) ; // overall_significance osig_l. add_set: add p-vals to multiset, overall sign. for graph
	//osig_r.add_set( pvals_r ) ;
	return ret ;
}

// evaluate one randomset
int* go_groups_hyper::calculate_rand( string &data, ostream *os ) 
{
	istringstream is( data.c_str() ) ;
//	int i = -1 ;
	int *ret = new int[10] ;
	for ( int i=0 ; i<10 ; ++i ) {
		ret[i] = 0 ;
	}

	vector<int> changed ;

	while ( is ) {
		int _changed ;
		is >> _changed ;
		changed.push_back( _changed ) ;
	}

	multiset<double> pvals_l, pvals_r ;
	// for every GO: get hyper p-vals
	for( unsigned int i = 0 ; i < names.size() ; ++i ) {
		// # genes annotated to root node
		double N = detected[root_idx] ;
		// # test genes annotated to root node 
		double M = changed[root_idx] ;
		// # genes annotated to current node
		double n = detected[i] ;
		// # test genes annotated to current node
		double x = changed[i] ;
		
		if ( n < cutoff ) continue ; 

		double prob_left;
		double prob_right ;

		// pval
		//				   ,- lower_tail
		// Steffi:
		prob_left  = R::phyper( x, M, N-M, n, 1, 0 ) ; 
		prob_right = R::phyper( x-1., M, N-M, n, 0, 0 ) ;

		pvals_l.insert( prob_left ) ;
		pvals_r.insert( prob_right ) ;

		//if ( os ) {
			//*os << names[i] << "\t" 
			   //<< N << "\t"
			   //<< n << "\t"
			   //<< M << "\t"
			   //<< x << "\t"
			   //<< prob_left << "\t"
			   //<< prob_right << endl ;
		//}
		// ret array contains number of significant GOs per p-cutoff
		if ( prob_left < 0.1 ) {
			ret[0]++ ;
			if ( prob_left < 0.05 ) {
				ret[1]++ ;
				if ( prob_left < 0.01 ) {
					ret[2]++ ;
					if ( prob_left < 0.001 ) {
						ret[3]++ ;
						if ( prob_left < 0.0001 ) {
							ret[4]++ ;
						}
					}
				}
			}
		}
		if ( prob_right < 0.1 ) {
			ret[5]++ ;
			if ( prob_right < 0.05 ) {
				ret[6]++ ;
				if ( prob_right < 0.01 ) {
					ret[7]++ ;
					if ( prob_right < 0.001 ) {
						ret[8]++ ;
						if ( prob_right < 0.0001 ) {
							ret[9]++ ;
						}
					}
				}
			}
		}
	}
	smallest_rand_p_l.insert( *(pvals_l.begin()) ) ;
	smallest_rand_p_r.insert( *(pvals_r.begin()) ) ;
	//osig_l.add_set( pvals_l ) ;
	//osig_r.add_set( pvals_r ) ;
	return ret ;
}

void go_groups_hyper::print_pvals( int nr_randsets, ostream &os ) {


	// loop over GOs and compute FWER
	for( unsigned int i = 0 ; i < names.size() ; ++i ) {
		if ( detected[i] >= cutoff ) { 
			// FWER: number of randsets with smallest p <= p-val 
			// for each group in dataset
			int n_l = 0 ; 
			multiset<double>::const_iterator it = smallest_rand_p_l.begin() ;
			while ( it != smallest_rand_p_l.end() && 
				*it <= data_pvals_l[i] + 1.0e-10 * data_pvals_l[i]) // NEW: add tolerance to account for float inaccuracy  
					n_l++, it++ ;
			int n_r = 0 ;
			it = smallest_rand_p_r.begin() ;
			while ( it != smallest_rand_p_r.end() && 
				*it <= data_pvals_r[i] + 1.0e-10 * data_pvals_r[i]) // NEW: add tolerance to account for float inaccuracy 
					n_r++, it++ ;
			// new: output higher precision for p-vals to check if FWER-order follows p-value-order
			os << std::setprecision(17) << names[i] << "\t"				
				<< data_pvals_l[i] << "\t" //p-values
				<< data_pvals_r[i] << "\t"
				<< static_cast<double>(n_l)/ static_cast<double>(nr_randsets) << "\t" // FWER_under 
				<< static_cast<double>(n_r)/ static_cast<double>(nr_randsets) << "\t" // FWER_over
				<< n_anno_expected[i] << "\t" << changed_data[i]  //NEW: exp and real number of candi genes
				<< endl;
		}
	} 
	
}
