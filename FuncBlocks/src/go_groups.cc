
#include "go_groups.h"
#include <cmath>
#include <map> 
#include <set> 
#include <fstream>

// steffi: Rcpp contains Rmath (eg. R::pnorm())
//#define MATHLIB_STANDALONE
//#include "../../include/Rmath.h"
#include <Rcpp.h>

using namespace std ;

go_groups::go_groups( string &groups, string &sgenes, int co, string root_go ) 
{
	cutoff = co ;

	istringstream is( groups.c_str() ) ;
	string name ;
	while ( is >> name ) {
		names.push_back( name ) ;
		// check.push_back( groups_to_check[name] ) ;
		if ( name == root_go ) root_idx = names.size() - 1 ; 
	}
	Rcpp::Rcerr << "GOs: " << names.size() << endl ;
	istringstream is2( sgenes.c_str() ) ;
	int nrgenes_per_group;
	while ( is2 >> nrgenes_per_group ) nr_of_genes.push_back( nrgenes_per_group ) ;

}


int* go_groups::calculate_data( string &data, double sum_nties, ostream *os ) 
{
	istringstream is( data.c_str() ) ;
	int *ret = new int[10] ;
	for ( int i=0 ; i<10 ; ++i ) {
		ret[i] = 0 ;
	}

	vector<double> ranksums;

	while ( is ) {
		double ranksum ;
		is >> ranksum ;
		ranksums.push_back( ranksum ) ;
	}


	data_pvals_l.resize( names.size() ) ;
	data_pvals_g.resize( names.size() ) ;

	multiset<double> pvals_less, pvals_greater ;

	for( unsigned int i = 0 ; i < names.size() ; ++i ) {
		
		data_pvals_l[i] = -1. ;
		data_pvals_g[i] = -1. ;

		double N = nr_of_genes[root_idx] ;
		double n = nr_of_genes[i] ;
		N -= n ;
		double ranksum = ranksums[i] ;
		
		double prob_greater ;
		double prob_less ;

		if ( n < cutoff ) continue ;

		double C = ranksum - ((n*(n+1.))/2.) ;
		double z = C - (n*N*0.5) ;
		double sigma = sqrt( (n*N/12.)*
				     ( (n+N+1.) - sum_nties / 
					((n+N)*(N+n-1.)) 
				     )
				   ) ; 
		double corr = -0.5 ;
		//steffi:
		prob_less = R::pnorm( (z-corr) / sigma, 0., 1., 1, 0 ) ;
		data_pvals_l[i] = prob_less ;
		corr = 0.5 ;
		//steffi:
		prob_greater = 1. - R::pnorm( (z-corr) / sigma, 0., 1., 1, 0 ) ;
		data_pvals_g[i] = prob_greater ;

		pvals_less.insert( prob_less ) ;
		pvals_greater.insert( prob_greater ) ;

		if ( os ) {
			*os << names[i] << "\t" 
			   << N << "\t"
			   << n << "\t"
		           << ranksum << endl ;
		}

		if ( prob_less < 0.1 ) {
			ret[0]++ ;
			if ( prob_less < 0.05 ) {
				ret[1]++ ;
				if ( prob_less < 0.01 ) {
					ret[2]++ ;
					if ( prob_less < 0.001 ) {
						ret[3]++ ;
						if ( prob_less < 0.0001 ) {
							ret[4]++ ;
						}
					}
				}
			}
		}
		if ( prob_greater < 0.1 ) {
			ret[5]++ ;
			if ( prob_greater < 0.05 ) {
				ret[6]++ ;
				if ( prob_greater < 0.01 ) {
					ret[7]++ ;
					if ( prob_greater < 0.001 ) {
						ret[8]++ ;
						if ( prob_greater < 0.0001 ) {
							ret[9]++ ;
						}
					}
				}
			}
		}
	}
	//less_sig.add_set( pvals_less ) ;
	//greater_sig.add_set( pvals_greater ) ;
	return ret ;
}
int* go_groups::calculate_rand( string &data, double sum_nties, ostream *os ) 
{
	istringstream is( data.c_str() ) ;
	int *ret = new int[10] ;
	for ( int i=0 ; i<10 ; ++i ) {
		ret[i] = 0 ;
	}

	vector<double> ranksums;

	while ( is ) {
		double ranksum ;
		is >> ranksum ;
		ranksums.push_back( ranksum ) ;
	}


	multiset<double> pvals_less, pvals_greater ;

	for( unsigned int i = 0 ; i < names.size() ; ++i ) {
		

		double N = nr_of_genes[root_idx] ;
		double n = nr_of_genes[i] ;
		N -= n ;
		double ranksum = ranksums[i] ;
		
		double prob_greater ;
		double prob_less ;

		if ( n < cutoff ) continue ;

		double C = ranksum - ((n*(n+1.))/2.) ;
		double z = C - (n*N*0.5) ;
		double sigma = sqrt( (n*N/12.)*
				     ( (n+N+1.) - sum_nties / 
					((n+N)*(N+n-1.)) 
				     )
				   ) ; 
		double corr = -0.5 ;
		prob_less = R::pnorm( (z-corr) / sigma, 0., 1., 1, 0 ) ;
//		if ( prob_less <= data_pvals_l[i] ) nr_rand_l[i]++ ;
		corr = 0.5 ;
		prob_greater = 1. - R::pnorm( (z-corr) / sigma, 0., 1., 1, 0 ) ;
//		if ( prob_greater <= data_pvals_g[i] ) nr_rand_g[i]++ ;

		pvals_less.insert( prob_less ) ;
		pvals_greater.insert( prob_greater ) ;

		if ( os ) {
			*os << names[i] << "\t" 
			   << prob_less << "\t"
			   << prob_greater << endl ;
		}

		if ( prob_less < 0.1 ) {
			ret[0]++ ;
			if ( prob_less < 0.05 ) {
				ret[1]++ ;
				if ( prob_less < 0.01 ) {
					ret[2]++ ;
					if ( prob_less < 0.001 ) {
						ret[3]++ ;
						if ( prob_less < 0.0001 ) {
							ret[4]++ ;
						}
					}
				}
			}
		}
		if ( prob_greater < 0.1 ) {
			ret[5]++ ;
			if ( prob_greater < 0.05 ) {
				ret[6]++ ;
				if ( prob_greater < 0.01 ) {
					ret[7]++ ;
					if ( prob_greater < 0.001 ) {
						ret[8]++ ;
						if ( prob_greater < 0.0001 ) {
							ret[9]++ ;
						}
					}
				}
			}
		}
	}
	smallest_rand_p_l.insert( *(pvals_less.begin()) ) ;
	//less_sig.add_set( pvals_less ) ;
	smallest_rand_p_g.insert( *(pvals_greater.begin()) ) ;
	//greater_sig.add_set( pvals_greater ) ;
	return ret ;
}


void go_groups::print_pvals( int nr_randsets, ostream &os ) {

	
	// vector<double> *fdr_qless = less_sig.fdr_qvals( 0 ) ; 
	// vector<double> *fdr_qgreater = greater_sig.fdr_qvals( 0 ) ; 

	//steffi: das unten war einkommentiert - fuer R-package keine FDR

	// map<double,double> *fdr_qless = less_sig.fdr_qvals( 0 ) ; 
	// map<double,double> *fdr_qgreater = greater_sig.fdr_qvals( 0 ) ; 

	for( unsigned int i = 0 ; i < names.size() ; ++i ) {
		if ( nr_of_genes[i] >= cutoff ) { 
			// FWER: number of randsets with smallest p <= p-val 
			// for each group in dataset
			int n_l = 0 ; 
			multiset<double>::const_iterator it = smallest_rand_p_l.begin() ;
			while ( it != smallest_rand_p_l.end() && 
				*it <= data_pvals_l[i] + 1.0e-10 * data_pvals_l[i]) // NEW: add tolerance to account for float inaccuracy
					n_l++, it++ ;
			int n_g = 0 ;
			it = smallest_rand_p_g.begin() ;
			while ( it != smallest_rand_p_g.end() && 
				*it <= data_pvals_g[i] + 1.0e-10 * data_pvals_g[i]) // NEW: add tolerance to account for float inaccuracy 
					n_g++, it++ ;
			os << names[i] << "\t" << data_pvals_l[i] << "\t"
				<< data_pvals_g[i] << "\t" 
				<< static_cast<double>(n_l)/
				   static_cast<double>(nr_randsets) << "\t" 
				<< static_cast<double>(n_g)/
				   static_cast<double>(nr_randsets) << endl;
				// steffi:   
				//<< (*fdr_qless)[data_pvals_l[i]] << "\t" 
				//<< (*fdr_qgreater)[data_pvals_g[i]] << endl ;
		}
	} 
	// steffi:
	//delete fdr_qless ;
	//delete fdr_qgreater ;

	// steffi: stoert im FUNC-package (out-Datei kein data.frame mehr)
	// steffi: falls man vllt. das einkommentiert braeuchte man vllt. noch ein delete fuer fdr_less und fdr_greater
	/*
	os << endl << endl 
           << "global-test-statistics (0 - 0.05): " << endl 
	   << less_sig.significance( 0, 0.05 ) << "\t" << greater_sig.significance( 0, 0.05 ) << endl ;
	os << endl ;
	

	vector<double> *fdr_less = less_sig.fdr( 0 ) ; 
	vector<double> *fdr_greater = greater_sig.fdr( 0 ) ; 

	os << "FDR" << endl ;
	os << "0.1\t0.05\t0.01\t0.001\t0.0001\t0.1\t0.05\t0.01\t0.001\t0.0001" <<  endl ;
	for ( int i = 0 ; i < 5 ; i++ ) os << (*fdr_less)[i] << "\t" ;
	for ( int i = 0 ; i < 5 ; i++ ) os << (*fdr_greater)[i] << "\t" ;
	os << endl ;
	os << endl ;
	*/

}
