
#include "go_groups_hyper.h"
#include <cmath>
#include <map> 
#include <set> 
#include <fstream>

// steffi: Rcpp contains Rmath (eg. R::pnorm())
//#define MATHLIB_STANDALONE
//#include "../../include/Rmath.h"
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
	while ( is >> name ) {
		is_det >> det ;
		is_ch >> ch ;
		//Rcpp::Rcerr << name << "\t" << det << "\t" << ch << endl ;
		names.push_back( name ) ;
		detected.push_back( det ) ;
		changed_data.push_back( ch ) ;
		if ( name == root_go ) root_idx = names.size() - 1 ; 
	}
	Rcpp::Rcerr << "GOs: " << names.size() << endl ;
//	Rcpp::Rcerr << detected.size() << endl ;
//	Rcpp::Rcerr << changed_data.size() << endl ;


}


int* go_groups_hyper::calculate_data( ostream *os ) 
{
//	int i = -1 ;
	int *ret = new int[10] ;
	for ( int i=0 ; i<10 ; ++i ) {
		ret[i] = 0 ;
	}

	vector<int> &changed = changed_data ;

	data_pvals_l.resize( names.size() ) ;
	data_pvals_r.resize( names.size() ) ;

	multiset<double> pvals_l, pvals_r ;

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
		//				   ,- lower_tail
		// steffi:
		prob_left  = R::phyper( x, M, N-M, n, 1, 0 ) ; 
		prob_right = R::phyper( x-1., M, N-M, n, 0, 0 ) ;
		//prob_left  = phyper( x, M, N-M, n, 1, 0 ) ; 
		//prob_right = phyper( x-1., M, N-M, n, 0, 0 ) ;
		data_pvals_l[i] = prob_left ;
		data_pvals_r[i] = prob_right ;

		pvals_l.insert( prob_left ) ;
		pvals_r.insert( prob_right ) ;

		if ( os ) {
			*os << names[i] << "\t" 
			   << N << "\t"
			   << n << "\t"
			   << M << "\t"
			   << x << "\t" << endl ;
		}

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
	osig_l.add_set( pvals_l ) ;
	osig_r.add_set( pvals_r ) ;
	return ret ;
}
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

	for( unsigned int i = 0 ; i < names.size() ; ++i ) {

		double N = detected[root_idx] ;
		double M = changed[root_idx] ;
		double n = detected[i] ;
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

		if ( os ) {
			*os << names[i] << "\t" 
			   << N << "\t"
			   << n << "\t"
			   << M << "\t"
			   << x << "\t"
			   << prob_left << "\t"
			   << prob_right << endl ;
		}

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
	osig_l.add_set( pvals_l ) ;
	osig_r.add_set( pvals_r ) ;
	return ret ;
}

void go_groups_hyper::print_pvals( int nr_randsets, ostream &os ) {

	// vector<double> *fdr_q_l = osig_l.fdr_qvals( 0 ) ; 
	// vector<double> *fdr_q_r = osig_r.fdr_qvals( 0 ) ; 

	/*steffi: fuer R-package keine FDR
	map<double,double> *fdr_q_l = osig_l.fdr_qvals( 0 ) ; 
	map<double,double> *fdr_q_r = osig_r.fdr_qvals( 0 ) ; 
	*/

	for( unsigned int i = 0 ; i < names.size() ; ++i ) {
		if ( detected[i] >= cutoff ) { 
			// FWER: number of randsets with smallest p <= p-val 
			// for each group in dataset
			int n_l = 0 ; 
			multiset<double>::const_iterator it = smallest_rand_p_l.begin() ;
			while ( it != smallest_rand_p_l.end() && 
				*it <= data_pvals_l[i] ) 
					n_l++, it++ ;
			int n_r = 0 ;
			it = smallest_rand_p_r.begin() ;
			while ( it != smallest_rand_p_r.end() && 
				*it <= data_pvals_r[i] ) 
					n_r++, it++ ;
			os << names[i] << "\t" << data_pvals_l[i] << "\t"
				<< data_pvals_r[i] << "\t" 
				<< static_cast<double>(n_l)/
				   static_cast<double>(nr_randsets) << "\t" 
				   << static_cast<double>(n_r)/ 
				   static_cast<double>(nr_randsets) << endl;
				/*steffi:
				<< (*fdr_q_l)[data_pvals_l[i]] << "\t" 
				<< (*fdr_q_r)[data_pvals_r[i]] << endl ;
				*/
		}
	} 
	/*steffi:
	delete fdr_q_l ;
	delete fdr_q_r ;
	*/
	// steffi: stoert im FUNC-package (out-Datei kein data.frame mehr)
	// steffi: falls man vllt. das einkommentiert braeuchte man vllt. noch ein delete fuer fdr_less und fdr_greater
	/*
	os << endl << endl 
           << "global-test-statistics (0 - 0.05): " << endl 
	    << osig_l.significance( 0, 0.05 ) << "\t" << osig_r.significance( 0, 0.05 ) << endl ;
	os << endl ;

//	ofstream ooo( "cdfs-left.txt" ) ;
//	osig_l.print_cdfs( ooo ) ;
//	ooo.close() ;

	vector<double> *fdr_l = osig_l.fdr( 0 ) ; 
	vector<double> *fdr_r = osig_r.fdr( 0 ) ; 

	os << "FDR" << endl ;
	os << "0.1\t0.05\t0.01\t0.001\t0.0001\t0.1\t0.05\t0.01\t0.001\t0.0001" <<  endl ;
	for ( int i = 0 ; i < 5 ; i++ ) os << (*fdr_l)[i] << "\t" ;
	for ( int i = 0 ; i < 5 ; i++ ) os << (*fdr_r)[i] << "\t" ;
	os << endl ;
	os << endl ;
	*/


}
