/*
 * FUNC - Functional Analysis of Gene Expression Data
 * Copyright (C) 2002  Bjoern Muetzel, Kay Pruefer
 * 
 * This program is modifiable/redistributable under the terms of the
 * GNU General Public License.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */


#include <time.h>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <set>
#include <vector>
#include <memory>

#include "go.h"
#include "go_graph_hyper.h"
#include "idmap.h"
#include "transitions.h"

#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;



//[[Rcpp::export]]
void hyper_randset(std::string detected, std::string changed ,int number_of_randomsets, std::string outfile, std::string term, std::string graph_path, std::string termtoterm, std::string root) 

{

// 1) Build GO-Graph using different files from go_date_termdb-tables.tar.gz
	
	// read term.txt
	std::ifstream terms( term.c_str() ) ;
	if ( ! terms ) {
		Rcerr << "Cannot open " << term << "." << endl ;
		//exit( 1 ) ;
	}
	idmap id_to_go( terms ) ;
	terms.close(  ) ;
	Rcerr << "Read " << id_to_go.size() << " terms." << endl ;

	// read graph_path.txt lesen
	std::ifstream transition_graph( graph_path.c_str() ) ;
	if ( ! transition_graph ) {
		Rcerr << "Cannot open " << graph_path << "." << endl ;
		//exit( 1 ) ;
	}
	string parent_go( root ) ;
	string parent_id = id_to_go.get_id_for_go( parent_go ) ;
	transitions trans( parent_id, transition_graph ) ;
	transition_graph.close(  ) ;
	Rcerr << "Found " << trans.size() << " nodes." << endl ;

	// read term2term lesen
	std::ifstream term2term( termtoterm.c_str() ) ;
	if ( ! term2term ) {
		Rcerr << "Cannot open " << termtoterm << "." << endl ;
		//exit( 1 ) ;
	}
	go_graph_hyper graph( trans, term2term, id_to_go ) ;
	term2term.close(  ) ;
	Rcerr << "Graph created." << endl ;
	
	// gos-object will be used to get one int* per GO and to print the results.
	go gos ;
	
	std::ostream *out ;
	out = new std::ofstream( outfile.c_str() ) ;

	std::ifstream in( detected.c_str() ) ;
	Rcerr << "Reading detectedfile... " 
			<< endl ;

// 2) Read expressed genes and their annotated GO-categories
	
	// gens == genes. This vector is a simple representation of the go tree.
	// every gene is 1 vector of int*, where the int* represents one go-node.
	vector< vector<int*> > gens;

	map<string,int> genename_to_index ;
	int index = 0 ;

	string line ;
	// in = detected = Allen:4005-file = one line per expr gene with genename and all its anno. GOs
	while ( in ) { 
		getline( in, line ) ;
		std::istringstream is( line.c_str() ) ; // take line as input stream
		string gen_name ;
		is >> gen_name ; // speichere erstes wort in line (genename) in gen_name
#ifdef DEBUG		
		Rcout << "gen_name: " << gen_name << endl ;
#endif
		vector<int*> gen_vec; // to store pointers to GOs for the current gene
		string go_name ;
		set<string> parents ;
		// go through all annotated GOs of an expressed gene
		while ( is >> go_name ) {
			// Get the names of all nodes that are parents of go_name
			graph.get_parents( go_name, &parents ) ;
		}
		for ( set<string>::const_iterator it = parents.begin() ; 
				it != parents.end() ; ++it ) 
		{
			// gos.add returns a unique int* for every string.
			gen_vec.push_back( gos.add( *it ) ) ; 
#ifdef DEBUG
			Rcout << "go: " << *it << endl ;
#endif
		}
		// if the gene is annotated, add it to the gens-vector (vector of vectors with GOs for all exp. genes)
		if ( gen_vec.size() > 0 ) {
			gens.push_back( gen_vec ) ;
			genename_to_index[gen_name] = index ;
			index++ ;	
		}	
	}
	Rcerr << "Found " << gens.size() << " usable entrys in " << detected
		<< " with " << gos.size() << " GOs" << endl ;
	
	*out << "Genes:\t" << gens.size() << endl ;
	*out << "GOs:\t" << gos.size() << endl ;
	
	// print line with names of GO-categories
	// print another line with sums of expressed genes annotated to the GOs (including background genes) 
	// also reset sums per GO to 0
	gos.print_names( *out ) ;

	
// 3) Read expressed test genes and add 1s to their annotated categories

	// read changedfile (infile data, one line per expressed TEST gene with genname and all GOs)
	// identical to in/detected but only for test genes
	// (Es wird nur der Genname gebraucht, der Rest der Zeile mit den GOs waere eigentlich nicht noetig)
	std::ifstream changed_in( changed.c_str() ) ;
	if ( ! changed_in ) {
		Rcerr << "Cannot open " << changed << endl ;
		//exit( 1 ) ;
	}
	int size_of_random_sets = 0 ; // number of expressed genes
	string line_s ;	
	while( changed_in ) {
		getline( changed_in, line_s ) ;
		std::istringstream is( line_s.c_str() ) ;
		string gen_name ;
		is >> gen_name ;
		Rcout << "annotated test gene: " << gen_name << endl;
		// genename_to_index: index of annotated GO-vec for a gene (gens[index])
		// if-statement should always be true because changed_in is subset of "in/detected" which was used to make genename_to_index 
		if ( genename_to_index.find( gen_name ) != genename_to_index.end() ) {
			// go through annotated GOs for the current gene and add 1
			for ( vector<int*>::iterator it = gens[genename_to_index[gen_name]].begin() ;
					it != gens[genename_to_index[gen_name]].end() ; ++it ) {
				(*(*it))++ ;
			}
			size_of_random_sets++ ;
		}
	}
	gos.print_sum( *out ) ;
	gos.clear() ;

	Rcerr << "Creating " << number_of_randomsets << " randomsets with "
			"size " << size_of_random_sets << endl ;


// 4) Create randomsets, randomly assign test genes among all annotated genes

	double max = gens.size() ; // gens.size = number of annotated genes
	Rcout << "Total number of annotated genes: " << max << endl;
	// forall randomsets
	for ( int i = 1 ; i <= number_of_randomsets ; ++i ) {
		// create randomset
		set<int> random_numbers ;		
		// size of random sets = anzahl annotierter test genes
		for ( int randi = 1 ; randi <= size_of_random_sets ; ++randi ) {
		
			int rand_num = static_cast<int>(max*R::runif(0,1)) ;
			// solange eine random number ziehen, bis eine kommt, die noch nicht im set ist
			while ( random_numbers.find( rand_num ) != random_numbers.end() ) {
				//Rcout << rand_num << " ";												
				rand_num = static_cast<int>(max*R::runif(0,1)) ;
			}
			
			//Rcout << rand_num << endl;
			
			random_numbers.insert( random_numbers.begin(), rand_num ) ;
		}		
		
		// reset all go-nodes, set all go_vec counters to zero
		gos.clear(  ) ;

		// add 1 to every GO that a randomly choosen gene is part of
		// go through all random numbers 
		for (set<int>::const_iterator it = random_numbers.begin() ;it != random_numbers.end() ; ++it){
			// gens == genes. This vector is a simple representation of the go tree.
			// every gene is 1 vector of int*, where the int represents one go-node.
			// gens[*it] = vector with all GOs that belong to the gene *it points at
			for ( vector<int*>::iterator it2 = gens[*it].begin() ;it2 != gens[*it].end() ; ++it2 ) {
				(*(*it2))++ ;
			}
		}
		random_numbers.clear() ;
		
		// print a line with the sums of annotaed genes for every go
		gos.print_sum( *out ) ;
	}
	delete out;
	
	Rcerr << "\rFinished" << endl ;
}
