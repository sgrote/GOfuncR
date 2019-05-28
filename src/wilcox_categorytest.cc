
/*
 * code modified from:
 * FUNC - Functional Analysis of Gene Expression Data
 * Copyright (C) 2002  Bjoern Muetzel, Kay Pruefer
 */

#include <vector>
#include <string>
#include <fstream>
#include "go_groups.h"

#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
void wilcox_category_test(std::string directory, int cut, std::string root, bool silent)
	
{

	/*************
         * parsing arguments, creating in and outstreams
         ************/
	string input = directory + "_randset_out";
	istream *in ;
	in = new ifstream( input.c_str() ) ;
	if ( ! *in ) {
		Rcpp::Rcerr << "Cannot open " << input << endl ;
	}

	string output = directory + "_category_test_out";
	ofstream out( output.c_str() ) ;
	if ( ! out ) {
		Rcpp::Rcerr << "Cannot open " << output << endl ;
	}
	
	string out_min_p = directory + "_min_p";
	ofstream min_p( out_min_p.c_str( )) ;
	if ( ! out ) {
		Rcpp::Rcerr << "Cannot open " << out_min_p << endl ;
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
	getline( *in, sites ) ; // sites (numbers of annotated genes per GO)

	if ( groups == "" || sites == "" ) {
		Rcpp::stop("Cant read Randomsets");
	}

	int co_genes_per_group = cut;

	/*************
         * go_groups handles parsing and analysis of dataset and randset lines
         ************/
	go_groups gos( groups, sites, co_genes_per_group, root ) ;

	string data ;
	getline( *in, data ) ; // real data (sums of scores of annotated genes)

	// returns number of significant groups for 0.1, 0.05, 0.01, 0.001, 0.0001
	int *realdata = gos.calculate_data( data, sum_nties);
	
	int sum_randdata[10] ;
	for ( int i=0 ; i < 10 ; ++i ) sum_randdata[i] = 0 ;
	int nr_groups_ge[10] ;
	for ( int i=0 ; i < 10 ; ++i ) nr_groups_ge[i] = 0 ;
	int num_randdata = 0 ;
	// randomsets
	if ( !silent ){
		Rcpp::Rcout << endl << "Evaluating randomsets: " << endl;
		Rcpp::Rcout << "No. of significant ontology nodes for" << endl;
		Rcpp::Rcout << "low ranks\t\t\t\thigh ranks" << endl ;
		Rcpp::Rcout << "of candidate genes at p-value thresholds" << endl ;
		Rcpp::Rcout << "0.1\t0.05\t0.01\t0.001\t0.0001\t0.1\t0.05\t0.01\t0.001\t0.0001" << endl ;
	}
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
		if ( !silent ){
			for ( int i=0 ; i < 10 ; ++i ) {
				Rcpp::Rcout << randdata[i] << "\t" ;
			}
			Rcpp::Rcout << "\n" ;
		}
		delete[] randdata ;
		num_randdata++ ;
	}
	// FWERs
	gos.print_pvals( num_randdata, out ) ;
	// save min_p to file for FWER-to-pval-interpolation in refinement
	gos.print_min_p( min_p ) ;
	
	// write summary to console
	if ( !silent ){
		Rcpp::Rcout << "Randomsets: " << num_randdata << endl << endl ;	
		// number of significant groups for 0.1, 0.05, 0.01, 0.001, 0.0001 for under- and over-rep -> 10 values
		Rcpp::Rcout << "Real data:" << endl ;
		Rcpp::Rcout << "No. of significant ontology nodes for" << endl;
		Rcpp::Rcout << "low ranks\t\t\t\thigh ranks" << endl ;
		Rcpp::Rcout << "of candidate genes at p-value thresholds" << endl ;
		Rcpp::Rcout << "0.1\t0.05\t0.01\t0.001\t0.0001\t0.1\t0.05\t0.01\t0.001\t0.0001" << endl ;
		for ( int i = 0 ; i < 10 ; ++i ) 
			Rcpp::Rcout << realdata[i] << "\t" ;
		Rcpp::Rcout << endl ;
		Rcpp::Rcout << "mean No. of significant groups in randomsets:" << endl ;
		for ( int i = 0 ; i < 10 ; ++i ) 
			Rcpp::Rcout << sum_randdata[i]/static_cast<double>(num_randdata) << "\t" ;
		Rcpp::Rcout << endl ;
		Rcpp::Rcout << "p value" << endl ;
		for ( int i = 0 ; i < 10 ; ++i ) 
			Rcpp::Rcout << nr_groups_ge[i]/static_cast<double>(num_randdata) << "\t" ;
		Rcpp::Rcout << endl << endl;
	}
	delete in;
	delete[] realdata;
}
