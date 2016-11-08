
/* vim: set ts=4 tw=80 cindent : */

#include "overall_sign.h" 
#include <fstream>
#include <algorithm>
#include <vector>

// used for iterating thru the pvals in the cdfs-matrix
#define MAX_PVALS static_cast<int>(1./precision)

overall_significance::overall_significance( double precision_ ) 
		: precision( precision_ ), nr_sets( 0 )  
{  }

void overall_significance::add_set( multiset<double> &pvals ) 
{

	cdfs_sets.push_back( pvals ) ; // vector<multiset<double> > cdfs_sets ;pvals: set of p-vals from original data or one randomset
	
	// vector to include into the cdfs
	vector<unsigned int> *vec = new vector<unsigned int>( MAX_PVALS+1, 0 ) ;
	
	// number of things with pvalue less than ... 
	int number_p=0 ;
	// idx of vector to add to the cdfs == pvalue 
	int idx=0;
	
	
	for ( multiset<double>::iterator it = pvals.begin() ; 
			it != pvals.end() ; ++it ) 
	{
		// skip gaps between pvalues and save 
		// number of pvals below for each p
		while( idx < index_for_pval(*it) ) {
			(*vec)[idx] = number_p ;
	   		idx++ ;
		}
		(*vec)[index_for_pval(*it)] = ++number_p ;
	}
	
	while( idx <= MAX_PVALS ) {
		(*vec)[idx] = number_p ;
		idx++ ;
	}
	
	cdfs.push_back( vec ) ;

	//steffi: vektor loeschen?
	delete vec ;
}

// debugging output
void overall_significance::print_cdfs( ostream &os ) 
{
	for ( int pval = 0 ; pval <= MAX_PVALS ; pval++ ) 
	{
		for ( int set = 0 ; set < cdfs.size() ; ++set )  
				os << (*cdfs[set])[pval] << "\t" ;
		os << endl ;
	}
	
}


// not used. average amount of things between 0 and 0.05
double overall_significance::alt_sign( int set_ ) 
{
	double average_data = 0.0 ;
	for ( int pval = 0 ; pval < MAX_PVALS*0.05 ; pval++ ) {
			int dataset = (*cdfs[set_])[pval] ;
			int nr_rand_le = 0 ;
			for ( int set = 0 ; set < cdfs.size() ; ++set ) {
					if ( set == set_ ) continue ;
					if ( (*cdfs[set])[pval] >= dataset ) 
							nr_rand_le++ ;
			}
			average_data += (  static_cast<double>(nr_rand_le)/
							   static_cast<double>(cdfs.size())  ) /
					        (MAX_PVALS*0.05) ;
	}
	return average_data ;
}

/* estimate the corrected p-value for each p and each set as the
   fraction of other datasets below. Search the maximum of the
   corrected pvalues between p=0 and co for each set. compare
   maximum of dataset number set_ to all other maxima and return
   fraction of other datasets that have a maxima at least as high 
   as the dataset */
   
double overall_significance::significance( int set_, double co ) 
{
	vector<int> maxima( cdfs.size(), -1 ) ; // maxima for each dataset
	vector<int> at_p( cdfs.size() ) ; // at which p-value

	// for all pvals
	for ( int pval = 0 ; pval < index_for_pval(co) ; pval++ ) {
		multimap<unsigned int, int> order ;
		// insert all values of sets @pval into multiset
		for ( int set = 0 ; set < cdfs.size() ; ++set ) { 
			order.insert( pair<unsigned int,int>( 
						(*(cdfs[set]))[pval], set ) ) ;
		}
		int sets_below=0 ;
		// it->first  = number of groups < pval
		// it->second = id of set
		// estimate for each how many sets have less groups at this pval 
		for ( multimap<unsigned int, int>::const_iterator it = 
			order.begin() ; it != order.end() ; ++it ) 
		{
			int tie = order.count( it->first ) ;
			int i = tie ;
			while ( i > 1 ) {
				i-- ;
				if ( maxima[it->second] < sets_below ) 
				{
					maxima[it->second] = sets_below ;
					at_p[it->second] = pval ;
				}
				it++ ;	
			}
			
			if ( maxima[it->second] < sets_below ) 
			{
				maxima[it->second] = sets_below ;
				at_p[it->second] = pval ;
			}
				
			sets_below+=tie ;	
		}
	}	

		
	// order maxima
	multimap<int, int> maxima_ordered ;
	for ( int i=0 ; i < maxima.size() ; i++ ) {
			maxima_ordered.insert( pair<int,int>( cdfs.size()-maxima[i], i ) ) ;
	}
	
	// estimate p-val
	int n = 0 ;
	for ( multimap<int, int>::const_iterator it = maxima_ordered.begin() ;
		it != maxima_ordered.end() ; ) 
	{
		int tie = maxima_ordered.count( it->first ) ;
		
		for ( int i = 0 ; i < tie ; i++ ) {
			if ( it->second == set_ ) {
				return ( static_cast<double>((tie+n)-1)/
			  		 	 static_cast<double>(cdfs.size()-1) ) ; 
			}
			it++ ;
		}
		n += tie ; 
	}
	// not reached
	return -1 ;
}

// return false discovery rate for cutoffs 0.0001, 0.001, 0.01, 0.05, 0.1
/* steffi: fdr nicht im R-packe
vector<double> *overall_significance::fdr( int set_ ) {

	if ( cdfs.size() == 0 ) return 0 ;
	// calculate fdrs for all p-values
	// if fdr < cutoff[i] add p to return-vector
	vector<double> *ret = new vector<double>( 5, -1. ) ;	// -1. -> no fdr found

	for ( int pval = 0 ; pval < MAX_PVALS ; pval++ ) {
		int sum_randsets = 0 ;	
		int dataset = 0 ;
		for ( int set = 0 ; set < cdfs.size() ; set++ ) {
			if ( set == set_ ) {
				dataset = (*cdfs[set])[pval] ;
				continue ;
			}
			sum_randsets += (*cdfs[set])[pval] ;
		}
		if ( dataset == 0 ) continue ;
		double V = static_cast<double>(sum_randsets) / 
				   static_cast<double>(cdfs.size() - 1 ) ;
		if ( V == 0. ) continue ;
		double FDR = V / static_cast<double>(dataset) ;
		if ( FDR < 0.1 ) {
			(*ret)[0] = pval*precision ;	
			if ( FDR < 0.05 ) {
				(*ret)[1] = pval*precision ;	
				if ( FDR < 0.01 ) {
					(*ret)[2] = pval*precision ;	
					if ( FDR < 0.001 ) {
						(*ret)[3] = pval*precision ;	
						if ( FDR < 0.0001 ) {
							(*ret)[4] = pval*precision ;	
						}	
					}	
				}
			}
		}
	}
	return ret ;
	
}

map<double,double> *overall_significance::fdr_qvals( int set_ ) 
{
	if ( cdfs_sets.size() == 0 ) return 0 ;
	map<double,double> *ret = new map<double,double> ;
	
	double m = (*cdfs[set_])[MAX_PVALS] ; // number of groups

	//cerr << "Number groups: " << m << endl ;
	
	// how many groups <= p for each randomset
	vector<int> nr_groups_per_randset( cdfs_sets.size(), 0 ) ;

	// iterators for all randomsets
	vector<multiset<double>::iterator> iterators ;
	for ( vector<multiset<double> >::const_iterator rand_it = cdfs_sets.begin() ;
					rand_it != cdfs_sets.end() ; ++rand_it ) 
	{
			iterators.push_back( rand_it->begin() ) ;
			
	}

	// go through all p-vals in dataset
	int r = 0 ; // number sign in dataset 
	for ( multiset<double>::const_iterator data_it = cdfs_sets[set_].begin() ;
					data_it != cdfs_sets[set_].end() ; ) 
	{
			double p = *data_it ;
			// count ties
			while( data_it != cdfs_sets[set_].end() && *data_it == p )
				r++, data_it++ ;
			
//			cerr << "Data: " << p << ", nr sign: " << r << endl ;
			
			// for all randomsets:
			double sum = 0. ;
			vector<int> rand_cumvals ;
			int i = 0 ;
			for ( vector<multiset<double> >::const_iterator rand_it = cdfs_sets.begin() ;
							rand_it != cdfs_sets.end() ; ++rand_it, ++i ) 
			{
					if ( i == set_ ) continue ; // ignore dataset 
//					cerr << "checking rand: " << i << endl ;
					multiset<double>::const_iterator &it = iterators[i] ;
					int x = 0 ;
					while ( it != rand_it->end() && *it <= p ) it++, x++ ;
//					cerr << "x = " << x << endl ;
					nr_groups_per_randset[i] += x ;
					sum += static_cast<double>(nr_groups_per_randset[i]) /
							( static_cast<double>(nr_groups_per_randset[i]) + r
							  			- m*p ) ;
					rand_cumvals.push_back( nr_groups_per_randset[i] ) ;
					//cerr << "... done: " << nr_groups_per_randset[i] << endl ;
			}
			nth_element( rand_cumvals.begin(),  
							rand_cumvals.begin()+int(rand_cumvals.size()*0.95), 
							rand_cumvals.end() ) ;
			double r_beta_star = *( rand_cumvals.begin() 
							        + int(rand_cumvals.size()*0.95) ) ;
			//cerr << "r_beta_star = " << r_beta_star << endl ;
			if ( r - r_beta_star >= p*m ) 
					(*ret)[p] = sum / static_cast<double>( cdfs_sets.size() - 1 ) ;
			else
					(*ret)[p] = -1. ;
			//cerr << "return: " << (*ret)[p] << endl ;
	}
				
	return ret ;
	
}
*/

/* steffi: das hier war schon auskommentiert
vector<double> *overall_significance::fdr_qvals( int set_ ) 
{
	if ( cdfs.size() == 0 ) return 0 ;
	vector<double> *ret = new vector<double>( MAX_PVALS+1, 1. ) ;
	
	double m = (*cdfs[set_])[MAX_PVALS] ; // number of groups

	for ( int pval = 0 ; pval < MAX_PVALS ; pval++ ) {
		double r = static_cast<double>((*cdfs[set_])[pval]) ; // number groups <= in dataset
		double p = pval*precision ; // pval between 0 and 1
		double sum = 0. ;
		vector<double> rand_pvals ;
		for ( int set = 0 ; set < cdfs.size() ; set++ ) {
			if ( set == set_ ) continue ;
			sum += static_cast<double>((*cdfs[set])[pval]) /
				( static_cast<double>((*cdfs[set])[pval]) + r - m*p ) ;
			rand_pvals.push_back( (*cdfs[set])[pval] ) ;
		}
		nth_element( rand_pvals.begin(),
			rand_pvals.begin()+int(rand_pvals.size()*0.95),
			rand_pvals.end() ) ;
		double r_beta_star = *( rand_pvals.begin()+ 
					int(rand_pvals.size()*0.95) ) ;
		if ( r - r_beta_star >= p*m ) 
			(*ret)[pval] = sum / static_cast<double>( cdfs.size()-1 ) ;
		else (*ret)[pval] = -1. ; // fwer should be used
			
	}
	return ret ;
}
*/
// index for p-vector 
int overall_significance::index_for_pval( double p ) 
{
	if ( p < 0. ) return 0 ;
	else if ( p > 1. ) return MAX_PVALS ;
	else return ( static_cast<int>(p/precision) ) ;
}


