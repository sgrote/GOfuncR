
#include "go_groups_conti.h"
#include <cmath>
#include <map> 
#include "fish_test.h"

#include <iomanip>      // std::setprecision
#include <Rcpp.h>

using namespace std ;

go_groups_conti::go_groups_conti( string &groups, istream *in, int co, string root ) 
{
	map<string, bool> groups_to_check ;
	if ( in ) {
		while ( *in ) {
			string go_gr ;
			*in >> go_gr ;
			int num_genes ;
			*in >> num_genes ;
			groups_to_check[go_gr] = (num_genes<co)?0:1 ;
		}	
	}
	istringstream is( groups.c_str() ) ;
	string name ;
	root_idx=-1 ;
	while ( is >> name ) {
		names.push_back( name ) ;
		check.push_back( groups_to_check[name] ) ;
		if ( name == root ) root_idx=names.size()-1 ;
	}
	if ( root_idx == -1 ) {
		Rcpp::Rcerr << "Cannot find GOID \"" << root << "\"." << endl ;
		Rcpp::stop("Failed to create go_groups_conti.\n") ;
	}
	//Rcpp::Rcout << "GOs: " << names.size() << endl ;
}


int* go_groups_conti::calculate_data( string &data, ostream *os ) 
{
	istringstream is( data.c_str() ) ;
	int *ret = new int[10] ;
	for ( int i=0 ; i<10 ; ++i ) {
		ret[i] = 0 ;
	}

	vector<int> ch_s ;
	vector<int> ch_ns ;
	vector<int> hh_s ;
	vector<int> hh_ns ;
	
	// p-values for (C/D > A/B) and (A/B > C/D)
	data_pvals_1.resize( names.size() ) ;
	data_pvals_2.resize( names.size() ) ;
	
	while ( is ) {
		int ch_s_, ch_ns_, hh_s_, hh_ns_ ;
		is >> ch_s_ ; // A  // synonymous and non-synonymous are switched in paper (Table 1)
		is >> ch_ns_ ; // B
		is >> hh_s_ ; // C
		is >> hh_ns_ ; // D
		ch_s.push_back( ch_s_ ) ;
		ch_ns.push_back( ch_ns_ ) ;
		hh_s.push_back( hh_s_ ) ;
		hh_ns.push_back( hh_ns_ ) ;
	}

	multiset<double> pvals_1, pvals_2 ;
	
	// loop over all nodes, get p-values	
	for ( int idx=0 ; idx < names.size() ; ++idx ) {

		data_pvals_1[idx] = -1. ;
		data_pvals_2[idx] = -1. ;

		if ( check[idx] == 0 ) continue ;

		double p_1 = 1. ;
		if ( hh_s[idx] == 0 && ch_s[idx] == 0 ) {
			if ( hh_ns[idx] < ch_ns[idx] ) p_1 = fisher_chi2( ch_s[idx], ch_ns[idx], hh_s[idx], hh_ns[idx] ) ;
			else p_1 = 1. ;
		} else if ( hh_s[idx] == 0 ) p_1 = 1. ;	
		else if ( ch_s[idx] == 0 ) p_1 = fisher_chi2( ch_s[idx], ch_ns[idx], hh_s[idx], hh_ns[idx] ) ;
		// if (D/C) >= (B/A), then p=1  (= A/B >= C/D)
		else if ( static_cast<double>(hh_ns[idx])/static_cast<double>(hh_s[idx]) 
                        >= static_cast<double>(ch_ns[idx])/static_cast<double>(ch_s[idx]) ) 
                   p_1 = 1. ;
		else 
                   p_1 = fisher_chi2( ch_s[idx], ch_ns[idx], hh_s[idx], hh_ns[idx] ) ;


		double p_2 = 1. ;
		if ( hh_s[idx] == 0 && ch_s[idx] == 0 ) {
			if ( hh_ns[idx] > ch_ns[idx] ) p_2 = fisher_chi2( ch_s[idx], ch_ns[idx], hh_s[idx], hh_ns[idx] ) ;
			else p_2 = 1. ;
		} else if ( ch_s[idx] == 0 ) p_2 = 1. ;	
		else if ( hh_s[idx] == 0 ) p_2 = fisher_chi2( ch_s[idx], ch_ns[idx], hh_s[idx], hh_ns[idx] ) ;
		// if (B/A) >= (D/C), then p=1 (= C/D >= A/B)
		else if ( static_cast<double>(ch_ns[idx])/static_cast<double>(ch_s[idx]) 
                        >= static_cast<double>(hh_ns[idx])/static_cast<double>(hh_s[idx]) ) 
                   p_2 = 1. ;
		else 
                   p_2 = fisher_chi2( ch_s[idx], ch_ns[idx], hh_s[idx], hh_ns[idx] ) ;

		data_pvals_1[idx] = p_1 ; // (C/D) > (A/B)
		data_pvals_2[idx] = p_2 ; // (A/B) > (C/D)

		pvals_1.insert( p_1 ) ;
		pvals_2.insert( p_2 ) ;

		//if ( os ) {
			//*os << names[idx] << '\t'
				//<< ch_s[idx] << '\t' 
				//<< ch_ns[idx] << '\t' 
				//<< hh_s[idx] << '\t'
				//<< hh_ns[idx] << endl ;
		//}
		// count number of nodes with p-value < 0.1, 0.05, ... 
		if ( p_1 < 0.1 ) {
			ret[0]++ ;
			if ( p_1 < 0.05 ) {
				ret[1]++ ;
				if ( p_1 < 0.01 ) {
					ret[2]++ ;
					if ( p_1 < 0.001 ) {
						ret[3]++ ;
						if ( p_1 < 0.0001 ) {
							ret[4]++ ;
						}
					}
				}
			}
		}
		if ( p_2 < 0.1 ) {
			ret[5]++ ;
			if ( p_2 < 0.05 ) {
				ret[6]++ ;
				if ( p_2 < 0.01 ) {
					ret[7]++ ;
					if ( p_2 < 0.001 ) {
						ret[8]++ ;
						if ( p_2 < 0.0001 ) {
							ret[9]++ ;
						}
					}
				}
			}
		}
	}
		

	//osig_1.add_set( pvals_1 ) ;
	//osig_2.add_set( pvals_2 ) ;

	return ret ;
}
int* go_groups_conti::calculate_rand( string &data, ostream *os ) 
{
	istringstream is( data.c_str() ) ;
	int *ret = new int[10] ;
	for ( int i=0 ; i<10 ; ++i ) {
		ret[i] = 0 ;
	}

	vector<int> ch_s ;
	vector<int> ch_ns ;
	vector<int> hh_s ;
	vector<int> hh_ns ;
	
	while ( is ) {
		int ch_s_, ch_ns_, hh_s_, hh_ns_ ;
		is >> ch_s_ ;
		is >> ch_ns_ ;
		is >> hh_s_ ;
		is >> hh_ns_ ;
		ch_s.push_back( ch_s_ ) ;
		ch_ns.push_back( ch_ns_ ) ;
		hh_s.push_back( hh_s_ ) ;
		hh_ns.push_back( hh_ns_ ) ;
	}

	multiset<double> pvals_1, pvals_2 ;
		
	for ( int idx=0 ; idx < names.size() ; ++idx ) {

		if ( check[idx] == 0 ) continue ;

		double p_1 = 1. ;
		if ( hh_s[idx] == 0 && ch_s[idx] == 0 ) {
			if ( hh_ns[idx] < ch_ns[idx] ) p_1 = fisher_chi2( ch_s[idx], ch_ns[idx], hh_s[idx], hh_ns[idx] ) ;
			else p_1 = 1. ;
		} else if ( hh_s[idx] == 0 ) p_1 = 1. ;	
		else if ( ch_s[idx] == 0 ) p_1 = fisher_chi2( ch_s[idx], ch_ns[idx], hh_s[idx], hh_ns[idx] ) ;
		else if ( static_cast<double>(hh_ns[idx])/static_cast<double>(hh_s[idx]) 
                        >= static_cast<double>(ch_ns[idx])/static_cast<double>(ch_s[idx]) ) 
                   p_1 = 1. ;
		else 
                   p_1 = fisher_chi2( ch_s[idx], ch_ns[idx], hh_s[idx], hh_ns[idx] ) ;


		double p_2 = 1. ;
		if ( hh_s[idx] == 0 && ch_s[idx] == 0 ) {
			if ( hh_ns[idx] > ch_ns[idx] ) p_2 = fisher_chi2( ch_s[idx], ch_ns[idx], hh_s[idx], hh_ns[idx] ) ;
			else p_2 = 1. ;
		} else if ( ch_s[idx] == 0 ) p_2 = 1. ;	
		else if ( hh_s[idx] == 0 ) p_2 = fisher_chi2( ch_s[idx], ch_ns[idx], hh_s[idx], hh_ns[idx] ) ;
		else if ( static_cast<double>(ch_ns[idx])/static_cast<double>(ch_s[idx]) 
                        >= static_cast<double>(hh_ns[idx])/static_cast<double>(hh_s[idx]) ) 
                   p_2 = 1. ;
		else 
                   p_2 = fisher_chi2( ch_s[idx], ch_ns[idx], hh_s[idx], hh_ns[idx] ) ;

		pvals_1.insert( p_1 ) ;
		pvals_2.insert( p_2 ) ;


		//if ( os ) {
			//*os << names[idx] << '\t'
				//<< ch_s[idx] << '\t' 
				//<< ch_ns[idx] << '\t' 
				//<< hh_s[idx] << '\t'
				//<< hh_ns[idx] << '\t'
				//<< p_1 << '\t'
				//<< p_2 << endl ;
		//}
		if ( p_1 < 0.1 ) {
			ret[0]++ ;
			if ( p_1 < 0.05 ) {
				ret[1]++ ;
				if ( p_1 < 0.01 ) {
					ret[2]++ ;
					if ( p_1 < 0.001 ) {
						ret[3]++ ;
						if ( p_1 < 0.0001 ) {
							ret[4]++ ;
						}
					}
				}
			}
		}
		if ( p_2 < 0.1 ) {
			ret[5]++ ;
			if ( p_2 < 0.05 ) {
				ret[6]++ ;
				if ( p_2 < 0.01 ) {
					ret[7]++ ;
					if ( p_2 < 0.001 ) {
						ret[8]++ ;
						if ( p_2 < 0.0001 ) {
							ret[9]++ ;
						}
					}
				}
			}
		}
	}
	// for family wise error rate
	smallest_rand_p_1.insert( *(pvals_1.begin()) ) ;
	smallest_rand_p_2.insert( *(pvals_2.begin()) ) ;
		
	//osig_1.add_set( pvals_1 ) ;
	//osig_2.add_set( pvals_2 ) ;

	return ret ;
}

void go_groups_conti::print_pvals( int nr_randsets, ostream &os ) {

	//// vector<double> *fdr_q_1 = osig_1.fdr_qvals( 0 ) ; 
	//// vector<double> *fdr_q_2 = osig_2.fdr_qvals( 0 ) ; 
	//map<double,double> *fdr_q_1 = osig_1.fdr_qvals( 0 ) ; 
	//map<double,double> *fdr_q_2 = osig_2.fdr_qvals( 0 ) ; 
	
	// loop over GOs and compute FWER
	for( unsigned int i = 0 ; i < names.size() ; ++i ) {
		if ( check[i] ) { 
			int n_1 = 0 ; 
			multiset<double>::const_iterator it = smallest_rand_p_1.begin() ;
			while ( it != smallest_rand_p_1.end() && 
				*it <= data_pvals_1[i] + 1.0e-10 * data_pvals_1[i])  // NEW: add tolerance
					n_1++,it++ ;
			int n_2 = 0 ;
			it = smallest_rand_p_2.begin() ;
			while ( it != smallest_rand_p_2.end() && 
				*it <= data_pvals_2[i] + 1.0e-10 * data_pvals_2[i])  // NEW: add tolerance
					n_2++,it++ ;
			os << std::setprecision(17)  << names[i] << "\t" 
				<< data_pvals_1[i] << "\t"
				<< data_pvals_2[i] << "\t" 
				<< static_cast<double>(n_1)/
				   static_cast<double>(nr_randsets) << "\t" 
				<< static_cast<double>(n_2)/
				   static_cast<double>(nr_randsets) << "\t"
				// TODO: hier noch was einfuegen so wie expected vs. real?  
				//<< (*fdr_q_1)[data_pvals_1[i]] << "\t" 
				//<< (*fdr_q_2)[data_pvals_2[i]] 
				<< endl ;
		}
	}

	//delete fdr_q_1 ;
	//delete fdr_q_2 ;

	//os << endl << endl 
           //<< "global-test-statistics (0 - 0.05): " << endl 
	   //<< osig_1.significance( 0, 0.05 ) << "\t" << osig_2.significance( 0, 0.05 ) << endl ;
	//os << endl ;

	//vector<double> *fdr_1 = osig_1.fdr( 0 ) ; 
	//vector<double> *fdr_2 = osig_2.fdr( 0 ) ; 

	//os << "FDR" << endl ;
	//os << "0.1\t0.05\t0.01\t0.001\t0.0001\t0.1\t0.05\t0.01\t0.001\t0.0001" <<  endl ;
	//for ( int i = 0 ; i < 5 ; i++ ) os << (*fdr_1)[i] << "\t" ;
	//for ( int i = 0 ; i < 5 ; i++ ) os << (*fdr_2)[i] << "\t" ;
	//os << endl ;
	//os << endl ;

}
