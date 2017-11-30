/*
 * code modified from:
 * FUNC - Functional Analysis of Gene Expression Data
 * Copyright (C) 2002  Bjoern Muetzel, Kay Pruefer
 */


#include <fstream>
#include <sstream>
#include <set>
#include <vector>

#include "idmap.h"
#include "transitions.h" 

#include "go_graph_conti.h"
#include "genes_conti.h"


#include <Rcpp.h>
using namespace Rcpp ;

//[[Rcpp::export]]
void conti_randset(std::string nodes_per_gene ,int number_of_randomsets, std::string directory, std::string root, bool silent)

{
	
	/*****************
         * read graph-structure and create graph
	 *******************/
	// read term.txt
	string term = directory + "_term.txt";
	std::ifstream terms( term.c_str() ) ;
	if ( ! terms ) {
		Rcpp::stop("Cannot open term.txt.\n"); 
	}
	idmap id_to_go( terms ) ;
	terms.close(  ) ;
	if ( !silent ){
		Rcpp::Rcout << "Read " << id_to_go.size() << " terms." << endl ;
	}
	// read graph_path
	string graph_path = directory + "_graph_path.txt";
	std::ifstream transition_graph( graph_path.c_str() ) ;
	if ( ! transition_graph ) {
		Rcpp::stop("Cannot open graph_path.txt.\n"); 
	}
	string parent_go( root ) ;
	string parent_id = id_to_go.get_id_for_go( parent_go ) ;
	transitions trans( parent_id, transition_graph ) ;
	transition_graph.close(  ) ;
	if ( !silent ){
		Rcpp::Rcout << "Found " << trans.size() << " nodes." << endl ;
	}
	// read term2term
	string termtoterm = directory + "_term2term.txt";
	std::ifstream term2term( termtoterm.c_str() ) ;
	if ( ! term2term ) {
		Rcpp::stop("Cannot open term2term.txt.\n"); 
	}
	go_graph_conti graph( trans, term2term, id_to_go ) ;
	term2term.close(  ) ;
	if ( !silent ){
		Rcpp::Rcout << "Graph created." << endl ;
	}
	
	

	/*****************
         * read gene information and annotate to graph
	 *******************/
	// read nodes_per_gene: one line per gene: gene | GO1 GO2 GO3
	std::ifstream annf( nodes_per_gene.c_str() ) ;
	if ( ! annf ) {
		Rcpp::stop("Cannot open nodes_per_gene.\n");
	}
	// read gene-scores: three columns: gene | count1 | count2
	string gene_scores = directory + "_infile-data";
	std::ifstream dataf( gene_scores.c_str() ) ;
	if ( ! dataf ) {
		Rcpp::stop("Cannot open infile-data.\n"); 
	}

	genes_conti gns( graph, annf, dataf ) ;
	
	if ( !silent ){
		Rcpp::Rcout << "Data and annotation file parsed." << endl ;
		Rcpp::Rcout << "Number of randomsets: " << number_of_randomsets << "." <<endl;	
		Rcpp::Rcout << "Computing randomsets..." << number_of_randomsets << "." <<endl;
	}
	

	/*****************
         * save #genes per go to extra file
	 *******************/
	// outfile for n_genes
	string outfile_ngene = directory + "_ngenes_per_go";
	ofstream out_n;
	out_n.open ( outfile_ngene.c_str() );
	graph.print_nr_genes( out_n ) ;
	out_n.close(  ) ;


	// outfile for randomsets
	string outfile = directory + "_randset_out";
	ofstream out;
	out.open ( outfile.c_str() );

	/*****************
         * save dataset to randomsetfile
	 *******************/
	graph.print_groups( out ) ;
	graph.print_vals( out ) ;

	/*****************
         * save random data to randomsetfile
	 *******************/
	for ( int i=1 ; i <= number_of_randomsets ; ++i ) {
		graph.clear_vals(  ) ;
		gns.create_random_set(  ) ;
		graph.print_vals( out ) ;
	}
	if ( !silent ){
		Rcpp::Rcout << "\rFinished" << endl ;
	}
	
}
