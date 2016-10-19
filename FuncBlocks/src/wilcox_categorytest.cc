
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include "go_groups.h"

// steffi:
#include <Rcpp.h>
using namespace Rcpp;

//steffi:
//int main( int argc, char *argv[] )
//[[Rcpp::export]]
void wilcox_category_test(std::string input, std::string output, int cut, std::string root)
	
{
	//steffi: veraltet: noch argument fuer profile-output: std::string output_profile, aber profile fuer R-package nicht noetig

	/*if ( argc != 6 ) {
		cerr << "Usage " << argv[0] << " randset outfile cutoff GO:ID profile" << endl ;
		//exit( 1 ) ;
	}*/
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
		in = new ifstream(  argv[1] ) ;
	}
	*/
	//steffi:
	if ( ! *in ) {
		Rcpp::Rcerr << "Cannot open " << input << endl ;
	}
	//steffi:
	ofstream out( output.c_str() ) ;
	if ( ! out ) {
		Rcpp::Rcerr << "Cannot open " << output << endl ;
	}

	// steffi:
	//ofstream *profile = new ofstream( output_profile.c_str() ) ;
	//if ( ! *profile ) profile = 0 ;

	string root_go ;
	{ 
		istringstream ppp( root.c_str() ) ;
		ppp >> root_go ;
	}

	/*************
         * start reading from randomset-file  
         ************/
	string s_sum_nties ;
	double sum_nties ;
	getline( *in, s_sum_nties ) ; // sum(NTIES^3 - NTIES)
	istringstream ntss( s_sum_nties.c_str() ) ;
	ntss >> sum_nties ;

	string groups, sites ;
	getline( *in, groups ) ; // GO IDs
	getline( *in, sites ) ; // sites

	if ( groups == "" || sites == "" ) {
		Rcpp::Rcerr << "Cant read Randomsets" << endl ;
		//exit( 1 ) ;
	}

	int co_genes_per_group = cut;
	// steffi: istringstream converts string to int, but already int handed in to R-function
	//int co_genes_per_group ;
	/*{
		istringstream cutoff( argv[3] ) ;
		cutoff >> co_genes_per_group ;
	}*/
	

	/*************
         * go_groups handles parsing and analysis of dataset and randset lines
         ************/
	go_groups gos( groups, sites, co_genes_per_group, root_go ) ;

	string data ;
	getline( *in, data ) ; // real data 

	// returns number of significant groups for 0.1, 0.05, 0.01, 0.001, 0.0001
	// steffi: generiert auch Profile-output -> outfile weglassen
	int *realdata = gos.calculate_data( data, sum_nties);
		//,profile ) ;
	
	//out << endl << endl ;

	int sum_randdata[10] ;
	for ( int i=0 ; i < 10 ; ++i ) sum_randdata[i] = 0 ;
	int nr_groups_ge[10] ;
	for ( int i=0 ; i < 10 ; ++i ) nr_groups_ge[i] = 0 ;
	int num_randdata = 0 ;
	// randomsets
	while ( *in ) {
		getline( *in, data ) ;
		if ( data == "" ) { break ; } 
		int *randdata = gos.calculate_rand( data, sum_nties ) ;
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
		//steffi
		delete[] randdata ;
		num_randdata++ ;
	}
	//gos.print_pvals( num_randdata, Rcpp::Rcout ) ;
	gos.print_pvals( num_randdata, out ) ;

	Rcpp::Rcout << "Randomsets: " << num_randdata << endl ;
	Rcpp::Rcout << "less\t\t\t\t\tgreater" << endl ;
	Rcpp::Rcout << "# sig. groups dataset" << endl ;
	for ( int i = 0 ; i < 10 ; ++i ) 
		Rcpp::Rcout << realdata[i] << "\t" ;
	Rcpp::Rcout << endl ;
	Rcpp::Rcout << "# sig. groups mean randomsets" << endl ;
	for ( int i = 0 ; i < 10 ; ++i ) 
		Rcpp::Rcout << sum_randdata[i]/static_cast<double>(num_randdata) << "\t" ;
	Rcpp::Rcout << endl ;
	Rcpp::Rcout << "# p value" << endl ;
	for ( int i = 0 ; i < 10 ; ++i ) 
		Rcpp::Rcout << nr_groups_ge[i]/static_cast<double>(num_randdata) << "\t" ;
	Rcpp::Rcout << endl ;
	
	// steffi: Fuer R-package nicht noetig:
	// write outfile
	/*
	out << "Randomsets: " << num_randdata << endl ;
	out << "less\t\t\t\t\tgreater" << endl ;
	out << "# sig. groups dataset" << endl ;
	for ( int i = 0 ; i < 10 ; ++i ) 
		out << realdata[i] << "\t" ;
	out << endl ;
	out << "# sig. groups mean randomsets" << endl ;
	for ( int i = 0 ; i < 10 ; ++i ) 
		out << sum_randdata[i]/static_cast<double>(num_randdata) << "\t" ;
	out << endl ;
	out << "# p value" << endl ;
	for ( int i = 0 ; i < 10 ; ++i ) 
		out << nr_groups_ge[i]/static_cast<double>(num_randdata) << "\t" ;
	out << endl ;
	*/

	// steffi: delete ifstream and realdata
	delete in;
	delete[] realdata;
}
