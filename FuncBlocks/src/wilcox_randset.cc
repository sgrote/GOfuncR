
#include <time.h>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <set>
#include <vector>
#include <memory>

#include "go_graph.h"
#include "idmap.h"
#include "transitions.h" 
#include "genes.h"

#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

#define MAX_LINE_LENGTH 20000

//[[Rcpp::export]]
void wilcox_randset(std::string nodes_per_gene ,int number_of_randomsets, std::string directory, std::string root) 
{
		
	/*****************
         * read graph-structure and create graph
	 *******************/
	// read term.txt
	string term = directory + "/term.txt";
	std::ifstream terms( term.c_str() ) ;
	if ( ! terms ) {
		Rcpp::stop("Cannot open term.txt.\n"); 
	}
	idmap id_to_go( terms ) ;
	terms.close(  ) ;
	Rcout << "Read " << id_to_go.size() << " terms." << endl ;
	
	// read graph_path	
	string graph_path = directory + "/graph_path.txt";
	std::ifstream transition_graph( graph_path.c_str() ) ;
	if ( ! transition_graph ) {
		Rcpp::stop("Cannot open graph_path.txt.\n"); 
	}
	string parent_go( root ) ;
	string parent_id = id_to_go.get_id_for_go( parent_go ) ;
	transitions trans( parent_id, transition_graph ) ;
	transition_graph.close(  ) ;
	Rcout << "Found " << trans.size() << " nodes." << endl ;
	
	// read term2term
	string termtoterm = directory + "/term2term.txt";
	std::ifstream term2term( termtoterm.c_str() ) ;
	if ( ! term2term ) {
		Rcpp::stop("Cannot open term2term.txt.\n"); 
	}
	go_graph graph( trans, term2term, id_to_go ) ;
	term2term.close(  ) ;
	Rcout << "Graph created." << endl ;

	/*****************
         * read gene information and annotate to graph
	 *******************/
	std::ifstream annf( nodes_per_gene.c_str() ) ;
	if ( ! annf ) {
		Rcpp::stop("Cannot open nodes_per_gene.\n");
	}
	string gene_scores = directory + "/infile-data";
	std::ifstream dataf( gene_scores.c_str() ) ;
	if ( ! dataf ) {
		Rcpp::stop("Cannot open infile-data.\n"); 
	}

	genes gns( graph, annf, dataf ) ;
	Rcout << "Data and annotation file parsed." << endl ;

	Rcout << "Number of randomsets: " << number_of_randomsets << "." <<endl;
	
	// steffi:
	Rcout << "Computing randomsets..." << number_of_randomsets << "." <<endl;

	string outfile = directory + "/randset_out";
	ofstream out;
	out.open ( outfile.c_str() );


	// force full length of all written floats...
	out.precision( 100 ) ; 

	/*****************
         * save original data, create and save randdata to randomsetfile, 
	 *******************/
	out << gns.sumnties() << endl ;
	graph.print_header( out ) ;
	graph.print_sumranks( out ) ;

	for ( int i=1 ; i <= number_of_randomsets ; ++i ) {
		graph.clear_genes(  ) ;
		gns.create_random_set(  ) ;
		graph.print_sumranks( out ) ;
	}

}
