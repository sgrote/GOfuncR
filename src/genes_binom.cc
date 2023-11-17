
#include "genes_binom.h"
#include <cstdlib>
#include <cstdio>
#include <Rcpp.h>
#include <random>
#include <algorithm>

genes_binom::genes_binom( go_graph_binom &graph, istream &annotation, istream &data ) 
{
	//srand( time(NULL) ) ;
	string line ;
	while ( annotation ) {
		getline( annotation, line ) ;
		istringstream is( line.c_str() ) ;
		string gene_name ;
		is >> gene_name ;

		if ( gene_name != "" ) {
			set<go_obj_binom*> parents ;
			string go_name ;
			while ( is >> go_name ) {
				graph.get_parents( go_name, &parents ) ;
			}
			if ( parents.size() > 0 ) {
				genemap[gene_name] = new gene_binom( gene_name, parents ) ;
			} else {
				Rcpp::Rcerr << gene_name << " not mapped.\n" ;
			}
		}
	}
	
	while( data ) {
		getline( data, line ) ;
		istringstream is( line.c_str() ) ;
		string gene_name ;
		is >> gene_name ;

		if ( genemap.find( gene_name ) != genemap.end() ) {
		
			double hka ; // first count
			is >> hka ;

			double cka ; // second count
			is >> cka ;

			if ( cka ==0. && hka == 0. ) {
				delete genemap[gene_name] ;
				genemap.erase( gene_name ) ;
			} else {
				genemap[gene_name]->add_hka_cka( static_cast<int>(hka), static_cast<int>(cka)) ;
				gene_vec.push_back( genemap[gene_name] ) ;
			}
		}
	}
}

genes_binom::~genes_binom(  ) 
{
	for ( vector<gene_binom*>::iterator it = gene_vec.begin() ; 
								it != gene_vec.end() ; ++it )  
	{
		delete *it ;
	}
}

void genes_binom::create_random_set(  ) 
{
  // std::random_device rd;
  std::mt19937 g(R::runif(0, 1e+05));
  std::shuffle( gene_vec.begin(), gene_vec.end(), g ) ;

	int i = 0 ;

	for ( map<string,gene_binom*>::const_iterator it = genemap.begin() ;
			it != genemap.end() ; ++it ) 
	{
		it->second->write_data_gos( gene_vec[i]->get_gos() ) ;
		i++ ;
	}
}
