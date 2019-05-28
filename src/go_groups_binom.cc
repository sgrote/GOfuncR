
#include "go_groups_binom.h"
#include <cmath>
#include <map> 
#include <fstream>

#include <iomanip>      // std::setprecision
#include <Rcpp.h>


using namespace std ;

// groups = names of groups, space seperated in string
// co = cutoff in number of genes
// in = number genes per go group
go_groups_binom::go_groups_binom( string &groups, istream *in, int co, string root_go ) 
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
	while ( is >> name ) {
		names.push_back( name ) ;
		check.push_back( groups_to_check[name] ) ;
		if ( name == root_go ) root_idx = names.size() - 1 ; 
	}
}

// data = left and right values whitespace separated for each group
int* go_groups_binom::calculate_data( string &data, ostream *os ) 
{
	// count number of nodes with p-value < 0.1, 0.05, ... 
	istringstream is( data.c_str() ) ;
	int i = -1 ;
	int *ret = new int[10] ;
	for ( int i=0 ; i<10 ; ++i ) {
		ret[i] = 0 ;
	}

	vector<int> human_ka ;
	vector<int> chimp_ka ;

	while ( is ) {
		i++ ;
		int c, h ;
		is >> c ;
		is >> h ;

		human_ka.push_back( h ) ;
		chimp_ka.push_back( c ) ;

	}

	// average for binomial test
	double p_binom = static_cast<double>(chimp_ka[root_idx])/
				static_cast<double>(chimp_ka[root_idx]+human_ka[root_idx]) ;

	// for overall_significance
	data_pvals_c.resize( names.size() ) ;
	data_pvals_h.resize( names.size() ) ;

	multiset<double> pvals_c, pvals_h ;
	// loop over all nodes, get p-value for over- and underrep
	for( unsigned int i = 0 ; i < names.size() ; ++i ) {

		data_pvals_c[i] = -1. ;
		data_pvals_h[i] = -1. ;

		int c = chimp_ka[i] ;
		int h = human_ka[i] ;

		if ( (c == 0) && (h == 0) ) continue ;
		if ( check[i] == 0 ) continue ;
		
		double pro_c = R::pbinom( c-1, c+h, p_binom, 0, 0 ) ;
		double pro_h = R::pbinom( h-1, c+h, (1.-p_binom), 0, 0 ) ;

		data_pvals_c[i] = pro_c ;
		data_pvals_h[i] = pro_h ;

		pvals_c.insert( pro_c ) ;
		pvals_h.insert( pro_h ) ;

		// count number of nodes with p-value < 0.1, 0.05, ... 
		if ( pro_c < 0.1 ) {
			ret[0]++ ;
			if ( pro_c < 0.05 ) {
				ret[1]++ ;
				if ( pro_c < 0.01 ) {
					ret[2]++ ;
					if ( pro_c < 0.001 ) {
						ret[3]++ ;
						if ( pro_c < 0.0001 ) {
							ret[4]++ ;
						}
					}
				}
			}
		}
		if ( pro_h < 0.1 ) {
			ret[5]++ ;
			if ( pro_h < 0.05 ) {
				ret[6]++ ;
				if ( pro_h < 0.01 ) {
					ret[7]++ ;
					if ( pro_h < 0.001 ) {
						ret[8]++ ;
						if ( pro_h < 0.0001 ) {
							ret[9]++ ;
						}
					}
				}
			}
		}
	}
	return ret ;
}

// evaluate one randomset
int* go_groups_binom::calculate_rand( string &data, ostream *os ) 
{
	istringstream is( data.c_str() ) ;
	int *ret = new int[10] ;
	for ( int i=0 ; i<10 ; ++i ) {
		ret[i] = 0 ;
	}

	vector<int> human_ka ;
	vector<int> chimp_ka ;

	while ( is ) {
		int c, h ;
		is >> c ;
		is >> h ;

		human_ka.push_back( h ) ;
		chimp_ka.push_back( c ) ;

	}

	// average for binomial test
	double p_binom = static_cast<double>(chimp_ka[root_idx])/
				static_cast<double>(chimp_ka[root_idx]+human_ka[root_idx]) ;

	multiset<double> pvals_c, pvals_h ;

	for( unsigned int i = 0 ; i < names.size() ; ++i ) {

		int c = chimp_ka[i] ;
		int h = human_ka[i] ;

		if ( (c == 0) && (h == 0) ) continue ;
		if ( check[i] == 0 ) continue ;
		
		// p-vals
		// P of getting at least c successes when trying c+h times with prob=(c/c+h in root)
		double pro_c = R::pbinom( c-1, c+h, p_binom, 0, 0 ) ; //P(C >= c)
		if ( ISNAN(pro_c) ) pro_c = 1. ; // happens if c < 0 
		double pro_h = R::pbinom( h-1, c+h, (1.-p_binom), 0, 0 ) ; //P(H >= h)
		if ( ISNAN(pro_h) ) pro_h = 1. ;

		pvals_c.insert( pro_c ) ;
		pvals_h.insert( pro_h ) ;

		// ret array contains number of significant GOs per p-cutoff
		if ( pro_c < 0.1 ) {
			ret[0]++ ;
			if ( pro_c < 0.05 ) {
				ret[1]++ ;
				if ( pro_c < 0.01 ) {
					ret[2]++ ;
					if ( pro_c < 0.001 ) {
						ret[3]++ ;
						if ( pro_c < 0.0001 ) {
							ret[4]++ ;
						}
					}
				}
			}
		}
		if ( pro_h < 0.1 ) {
			ret[5]++ ;
			if ( pro_h < 0.05 ) {
				ret[6]++ ;
				if ( pro_h < 0.01 ) {
					ret[7]++ ;
					if ( pro_h < 0.001 ) {
						ret[8]++ ;
						if ( pro_h < 0.0001 ) {
							ret[9]++ ;
						}
					}
				}
			}
		}
	}
	
	// for family wise error rate
	smallest_rand_p_c.insert( *(pvals_c.begin()) ) ;
	smallest_rand_p_h.insert( *(pvals_h.begin()) ) ;

	return ret ;
}


void go_groups_binom::print_pvals( int nr_randsets, ostream &os ) {
	
	// loop over GOs and compute FWER
	for( unsigned int i = 0 ; i < names.size() ; ++i ) {
		if ( check[i] == 1 ) { 
			int n_c = 0 ; 
			multiset<double>::const_iterator it = smallest_rand_p_c.begin() ;
			while ( it != smallest_rand_p_c.end() && 
				*it <= data_pvals_c[i] + 1.0e-10 * data_pvals_c[i]) // add tolerance to account for float inaccuracy  
					n_c++, it++ ;
			int n_h = 0 ;
			it = smallest_rand_p_h.begin() ;
			while ( it != smallest_rand_p_h.end() && 
				*it <= data_pvals_h[i] + 1.0e-10 * data_pvals_h[i]) // add tolerance to account for float inaccuracy  
					n_h++, it++ ;
			os << std::setprecision(17) << names[i] << "\t" 				
				<< data_pvals_c[i] << "\t" //p-values
				<< data_pvals_h[i] << "\t"
				<< static_cast<double>(n_c)/ static_cast<double>(nr_randsets) << "\t"  //FWER high A
				<< static_cast<double>(n_h)/ static_cast<double>(nr_randsets) << "\t"  //FWER high B
				<< endl ;
		}
	} 

}

void go_groups_binom::print_min_p( ostream &os ) {
	
	multiset<double>::const_iterator it_l = smallest_rand_p_c.begin() ;
	multiset<double>::const_iterator it_r = smallest_rand_p_h.begin() ;
	while ( it_l != smallest_rand_p_c.end() ) { 
		os << std::setprecision(17) << *it_l << "\t" << *it_r << endl;
		it_l++, it_r++ ;
	}
}
