//steffi:
//using namespace std ;

#include <time.h>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <set>
#include <vector>
#include <memory>

// steffi:
#include "go_graph.h"
#include "idmap.h"
#include "transitions.h" 
#include "genes.h"

#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

#define MAX_LINE_LENGTH 20000

// steffi:
//int main( int argc, char *argv[] ) 

//[[Rcpp::export]]
void wilcox_randset(std::string nodes_per_gene, std::string gene_scores ,int number_of_randomsets, std::string outfile, std::string term, std::string graph_path, std::string termtoterm, std::string root) 
{
	/*steffi:
		if (  argc != 9 ) {
		cerr << "Usage: " << argv[0] << " infile_ann "
				"infile_data number_of_sets outfile term.txt" << endl
				<< "       graph_path.txt term2term.txt GO_ID" << endl ;
		//exit( 1 ) ;
	}*/
		
	/*****************
         * read graph-structure and create graph
	 *******************/
	// steffi: term.txt lesen
	std::ifstream terms( term.c_str() ) ;
	if ( ! terms ) {
		Rcerr << "Cannot open " << term << "." << endl ;
		//exit( 1 ) ;
	}
	idmap id_to_go( terms ) ;
	terms.close(  ) ;
	Rcerr << "Read " << id_to_go.size() << " terms." << endl ;
	
	// steffi: graph_path.txt lesen
	std::ifstream transition_graph( graph_path.c_str() ) ;
	if ( ! transition_graph ) {
		Rcerr << "Cannot open " << graph_path << "." << endl ;
		//exit( 1 ) ;
	}
	//steffi:
	string parent_go( root ) ;
	string parent_id = id_to_go.get_id_for_go( parent_go ) ;
	transitions trans( parent_id, transition_graph ) ;
	transition_graph.close(  ) ;
	Rcerr << "Found " << trans.size() << " nodes." << endl ;
	
	// steffi: term2term lesen
	std::ifstream term2term( termtoterm.c_str() ) ;
	if ( ! term2term ) {
		Rcerr << "Cannot open " << termtoterm << "." << endl ;
		//exit( 1 ) ;
	}
	go_graph graph( trans, term2term, id_to_go ) ;
	term2term.close(  ) ;
	Rcerr << "Graph created." << endl ;

	/*****************
         * read gene information and annotate to graph
	 *******************/
    //steffi: 
	std::ifstream annf( nodes_per_gene.c_str() ) ;
	if ( ! annf ) {
		Rcerr << "Cannot open " << nodes_per_gene << "." << endl ;
		//exit( 1 ) ;
	}
	//steffi:
	std::ifstream dataf( gene_scores.c_str() ) ;
	if ( ! dataf ) {
		Rcerr << "Cannot open " << gene_scores << "." << endl ;
		//exit( 1 ) ;
	}

	genes gns( graph, annf, dataf ) ;
	Rcerr << "Data and annotation file parsed." << endl ;
	// steffi:
	//int number_of_randomsets ;
	//{
		//istringstream is( argv[3] ) ;
		//is >> number_of_randomsets ;
	//}
	Rcout << "Number of randomsets: " << number_of_randomsets << "." <<endl;
	
	// steffi:
	Rcout << "Computing randomsets..." << number_of_randomsets << "." <<endl;

	// steffi: Udos Vorschlag fuer output
	// Original:
	/*
	ostream *out ;
	if ( outfile == "-" ) {
	/ $cout: address of cout	
		out = &cout ;
	} else {
		out = new ofstream( outfile.c_str() ) ;
		if ( !*out ) {
			Rcerr << "Cannot open output file " << outfile << endl ;
			////exit( 1 ) ;
		}
	}
	
	out->precision( 100 ) ; 
	
	*out << gns.sumnties() << endl ;
	graph.print_header( *out ) ;
	graph.print_sumranks( *out ) ;

	for ( int i=1 ; i <= number_of_randomsets ; ++i ) {
		graph.clear_genes(  ) ;
		gns.create_random_set(  ) ;
		graph.print_sumranks( *out ) ;		
	}
	*/
	
	ofstream out;
	//if (outfile == "-"){
	//	out.ios::rdbuf( cout.rdbuf() );
	//} else {
	  	out.open ( outfile.c_str() );
	//}	
	/* steffi: Udo:Jetzt ist allerdings 'out' kein Zeiger mehr, sondern ein Objekt.
	Folglich wird danach aus jedem '*out' ein gewoehnliches 'out', aus jedem
	'out->' ein 'out.' und aus jedem 'out' ein '&out'.*/

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
