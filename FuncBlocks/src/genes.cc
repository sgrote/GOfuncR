
#include "genes.h"
//#include <time.h>
#include <cstdlib>
#include <cstdio>
//#include <cmath>
//#define MAX_LINE_LENGTH 10000
#include <Rcpp.h>

genes::genes( go_graph &graph, istream &annotation, istream &data ) 
{
//	srand( time(NULL) ) ;
//	srand( 100 ) ;
	// annotation = root-file: one line per gene: gene | GO1 GO2 GO3
	string line ;
	while ( annotation ) {
		getline( annotation, line ) ;
		istringstream is( line.c_str() ) ;
		string gene_name ;
		is >> gene_name ;

		if ( gene_name != "" ) {
			set<go_obj*> parents ;
			string go_name ;
			while ( is >> go_name ) {
				graph.get_parents( go_name, &parents ) ;
			}
			if ( parents.size() > 0 ) {
				// "gene" has gene-name and annotated-GOs including parent GOs
				genemap[gene_name] = new gene( gene_name, parents ) ; 
			} else {
				Rcpp::Rcerr << gene_name << " not mapped.\n" ;
			}
		}
	}
	Rcpp::Rcerr << "Annotated " << genemap.size() << " genes." << endl ;
	
	multimap<double, gene*> genes_ranked ;
	// data: two columns: gene | score
	while( data ) {
		getline( data, line ) ;
		istringstream is( line.c_str() ) ;
		string gene_name ;
		is >> gene_name ;

		if ( genemap.find( gene_name ) != genemap.end() ) {
		
			double data ; // score
			is >> data ;

			gene_vec.push_back( genemap[gene_name] ) ;
			// score and gene (with anno GOs)
			genes_ranked.insert(pair<double, gene*>(data,genemap[gene_name]) ) ;
		}
	}
	// rank
	int i = 1 ;
	sum_nties = 0. ;
	// loop over unique scores, set ranks for every gene in "genes_ranked" (gene* = gene and GOs)
	for ( multimap< double, gene* >::iterator it = genes_ranked.begin() ; it != genes_ranked.end() ;  ){ 
		int equal = genes_ranked.count( it->first ) ; // how many have the same score?
		//Rcpp::Rcout << "\nequal:\n" << equal << "\nscore:\n" << it->first << std::endl;
		if ( equal > 1 ) {
			sum_nties += pow(static_cast<double>(equal), 3.) // equal^3 - equal
					-static_cast<double>(equal) ;
			// for every gene with the same score 
			// set the same averaged rank (scores: 2,4,4,6 -> ranks: 1,2.5,2.5,4)	
			double rank = static_cast<double>(i) + static_cast<double>(equal-1)/2. ;
			for ( int i2 = 0 ; i2 < equal ; i2++ ) {
				//Rcpp::Rcout << "tied rank:\n" << rank << std::endl;
				it->second->set_rank( rank ) ; //gene.set_rank
				i++ ;
				++it ;
			}
		} else {
			// un-tied: rank = position of gene in the set
			//Rcpp::Rcout << "normal rank:\n" << i << std::endl;
			it->second->set_rank( i ) ;
			i++ ;
			++it ;
		}
	}
	// genemap = map<gene-name, gene*>; class gene = (name, set<go_obj*>)
	for ( map<string,gene*>::const_iterator it = genemap.begin() ; it != genemap.end() ; ++it ) {
		it->second->write_to_gos( ) ; 
		// gene.write_to_gos -> go-obj.add_gene(this-gene) -> genes.push_back(gene) => vector<gene*> genes
	}
}

genes::~genes(  ) 
{
	for ( map<string,gene*>::const_iterator it = genemap.begin() ; it != genemap.end() ; ++it )	{
		delete it->second ;
	}
}

#include <algorithm> 

class c_prng {
	public:
// steffi:
//		int operator()( int n ) { return( rand() % n ); } ;
		int operator()( int n ) { return( R::runif(0,1) * n ); } ;
} ;

void genes::create_random_set(  ) 
{
	c_prng rng ;
        random_shuffle( gene_vec.begin(), gene_vec.end(), rng ) ;

        int i = 0 ;

        for ( map<string,gene*>::const_iterator it = genemap.begin() ;
                        it != genemap.end() ; ++it )
        {
                it->second->write_to_gos( gene_vec[i]->get_gos() ) ;
                i++ ;
        }
}
