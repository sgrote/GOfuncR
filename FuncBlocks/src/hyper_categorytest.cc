
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include "go_groups_hyper.h"

// steffi:
#include <Rcpp.h>
using namespace Rcpp;

//steffi:
//int main( int argc, char *argv[] )
//[[Rcpp::export]]
void hyper_category_test(std::string input, std::string output, int cutoff, std::string root)

{
	/*steffi: veraltet: noch argument fuer profile-output: std::string output_profile, aber profile fuer R-package nicht noetig

	if ( argc != 6 ) {
		cerr << "Usage " << argv[0] << " randset outfile profile cutoff GO:ID" << endl ;
		exit( 1 ) ;
	}
	*/
	/*************
         * parsing arguments, creating in and outstreams
         ************/
	// steffi: in R-package wird nur Datei gelesen, und so kann einfach ein "delete" benutzt werden
	istream *in ;
	in = new ifstream( input.c_str() ) ;
	/*
	if ( input == "-" ) {
		in = &cin ;
	} else {
		in = new ifstream( argv[1] ) ;
	}
	*/
	//steffi:
	if ( ! *in ) {
		Rcpp::Rcerr << "Cannot open " << input << endl ;
	}
	//steffi:
	ofstream out( output.c_str( )) ;
	if ( ! out ) {
		Rcpp::Rcerr << "Cannot open " << output << endl ;
	}

	/* steffi: istringstream converts string to int, but already int handed in to R-function
	int cutoff=1 ;
	{
		istringstream cutoff_s( argv[4] ) ;
		cutoff_s >> cutoff ;
	}
	*/

	string root_go ;
	{	
		//steffi:
		//istringstream ppp( argv[5] ) ;
		istringstream ppp( root.c_str() ) ;
		ppp >> root_go ;
	}


	/*************
         * start reading from randomset-file  
         ************/
	// ignore header-lines
	string dummy ;
	getline( *in, dummy ) ;  // #Genes
	getline( *in, dummy ) ;  // #annotated GOs

	string groups ;	
	getline( *in, groups ) ; // GO IDs
	if ( groups == "" ) {
		//steffi:
		Rcpp::Rcerr << "Error reading randomsets" << endl ;
		//exit( 1 ) ;
	}
	
	string detected, changed ;
	// ignore detected-values (reading it from profile instead)
	getline( *in, detected ) ;  // #annotated genes (candidate and background)
	getline( *in, changed ) ;  // #annotated candidate genes

	/*************
         * go_groups handles parsing and analysis of dataset and randset lines
         ************/
	go_groups_hyper gos( groups, detected, changed, root_go, cutoff ) ;

	// returns number of significant groups for 0.1, 0.05, 0.01, 0.001, 0.0001
	// for under and over respresentation -> 10 values, also printed to console below
	int *realdata = gos.calculate_data( ) ;

//	int *realdata = gos.calculate_data( &profile ) ;
//	out << endl << endl ;
	

	int sum_randdata[10] ; // sum significant groups across all randomsets 
	for ( int i=0 ; i < 10 ; ++i ) sum_randdata[i] = 0 ;
	int nr_groups_ge[10] ; // sum of randsets where randdata[i] >= realdata[i] 
	for ( int i=0 ; i < 10 ; ++i ) nr_groups_ge[i] = 0 ;
	int num_randdata = 0 ;
	// randomsets
	Rcpp::Rcout << endl << "Evaluating randomsets: " << endl;
	Rcpp::Rcout << "No. of significant ontology nodes for" << endl;
	Rcpp::Rcout << "under-representation\t\t\tover-representation" << endl ;
	Rcpp::Rcout << "of candidate genes at p-value thresholds" << endl ;
	Rcpp::Rcout << "0.1\t0.05\t0.01\t0.001\t0.0001\t0.1\t0.05\t0.01\t0.001\t0.0001" << endl ;
	string data ;
	while ( *in ) {
		getline( *in, data ) ;
		if ( data == "" ) { break ; } 
		int *randdata = gos.calculate_rand( data ) ;
		for ( int i=0 ; i<10 ; ++i ) {
			sum_randdata[i] += randdata[i] ;
			if ( randdata[i] >= realdata[i] ) {
				nr_groups_ge[i]++ ;
			}
		}
		for ( int i=0 ; i < 10 ; ++i ) {
			Rcpp::Rcout << randdata[i] << "\t" ;
		}
		Rcpp::Rcout << "\n" ;
		//steffi:
		delete[] randdata ;
		num_randdata++ ;
	}
	// gos.print_pvals( num_randdata, Rcpp::Rcout ) ;
	gos.print_pvals( num_randdata, out ) ;

	// write summary to console
	Rcpp::Rcout << "Randomsets: " << num_randdata << endl << endl ;	
	// number of significant groups for 0.1, 0.05, 0.01, 0.001, 0.0001 for under- and over-rep -> 10 values
	Rcpp::Rcout << "Real data:" << endl ;
    Rcpp::Rcout << "No. of significant ontology nodes for" << endl;
	Rcpp::Rcout << "under-representation\t\t\tover-representation" << endl ;
	Rcpp::Rcout << "of candidate genes at p-value thresholds" << endl ;
	Rcpp::Rcout << "0.1\t0.05\t0.01\t0.001\t0.0001\t0.1\t0.05\t0.01\t0.001\t0.0001" << endl ;
	for ( int i = 0 ; i < 10 ; ++i ) 
		Rcpp::Rcout << realdata[i] << "\t" ;
	Rcpp::Rcout << endl ;
	// Mean number of significant nodes across all randomsets for different p-value cutoffs 0.1, 0.05, 0.01, 0.001, 0.0001; under- and overrep
	//Rcpp::Rcout << "# sig. groups mean randomsets" << endl ;
	Rcpp::Rcout << endl << "mean No. of significant groups in randomsets:" << endl ;

	for ( int i = 0 ; i < 10 ; ++i ) 
		Rcpp::Rcout << sum_randdata[i]/static_cast<double>(num_randdata) << "\t" ;
	Rcpp::Rcout << endl ;
	// proportion of randomsets where number of significant nodes at p 0.1, 0.05, 0.01, 0.001, 0.0001 is >= in original data
	Rcpp::Rcout << "p value" << endl ;
	for ( int i = 0 ; i < 10 ; ++i ) 
		Rcpp::Rcout << nr_groups_ge[i]/static_cast<double>(num_randdata) << "\t" ;
	Rcpp::Rcout << endl << endl ;


	// steffi: delete ifstream and realdata
	delete in;
	delete[] realdata;
}
