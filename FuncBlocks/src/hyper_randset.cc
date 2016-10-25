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


//#include <time.h>
//#include <cstdlib>
//#include <cstdio>
#include <fstream>
#include <sstream>
#include <set>
#include <vector>
#include <memory>
#include <map> 		
#include <iostream> 

#include "go.h"
#include "go_graph_hyper.h"
#include "idmap.h"
#include "transitions.h"

#include "structures.h"
#include "read_bed.h"
#include "blocks.h"
#include "roll.h"
#include "ran_genelen.h"

#include <Rcpp.h>


using namespace Rcpp;



//[[Rcpp::export]]
void hyper_randset(std::string all_genes, int number_of_randomsets, std::string directory, std::string root, std::string mod){
	
// 1) Build GO-Graph using different files from go_date_termdb-tables.tar.gz
	
	// read term.txt
	string term = directory + "/term.txt";
	std::ifstream terms( term.c_str() ) ;
	if ( ! terms ) {
		Rcpp::stop("Cannot open term.txt.\n"); 
	}
	idmap id_to_go( terms ) ;
	terms.close(  ) ;
	Rcout << "Read " << id_to_go.size() << " terms." << endl ;

	// read graph_path.txt
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

	// read term2term.txt
	string termtoterm = directory + "/term2term.txt";
	std::ifstream term2term( termtoterm.c_str() ) ;
	if ( ! term2term ) {
		Rcpp::stop("Cannot open term2term.txt.\n");
	}
	go_graph_hyper graph( trans, term2term, id_to_go ) ;
	term2term.close(  ) ;
	Rcout << "Graph created." << endl ;
	


// 2) Read expressed genes, their positions and their annotated GO-categories
		// gos-object will be used to get one int* per GO and to print the results.
	go gos ;	

	// maybe rename detected, changed to all_genes and cadidate_genes oder so	
	std::ifstream in( all_genes.c_str() ) ;
	Rcout << "Reading all_genes file... " << endl ;
			
	// gens == genes. This vector is a simple representation of the go tree.
	// every gene is 1 vector of int*, where the int* represents one go-node.
	vector< vector<int*> > gens;
	// vector of structures to store coordinates of genes from Allen:4005
	vector<gen_pos_str> genes_pos;
	// genename to index for gens
	map<string,int> genename_to_index ;
	int index = 0 ;
	long total_length = 0; // for gene_len option	
	// one line per expr gene: gene | (chrom) | (start) | (end) | GO1 GO2 GO3 
	string line ;	
	while ( getline(in, line) ) { 
		std::istringstream is( line.c_str() ) ; 
		// store name and optionally position of current gene
		gen_pos_str gen_pos; 
		is >> gen_pos.name;
		if(mod != "classic"){	
			is >> gen_pos.chrom >> gen_pos.start >> gen_pos.end;		
		}
		if(mod == "gene_len"){
			total_length += gen_pos.end - gen_pos.start;
			gen_pos.cumu_len = total_length;
		}
		genes_pos.push_back(gen_pos);
		// go through all annotated GOs of an expressed gene and store pointers to GOs
		vector<int*> gen_vec; 
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
			//Rcout << "go: " << *it << endl ;

		}
		// if the gene is annotated, add it to the gens-vector (vector of vectors with GOs for all exp. genes)
		if ( gen_vec.size() > 0 ) {
			gens.push_back( gen_vec ) ;
			genename_to_index[gen_pos.name] = index ;
			index++ ;	
		}	
	}
	in.close();
	
	int n_background = gens.size();
	
	Rcout << "Found " << n_background << " usable entrys in " << all_genes
		<< " with " << gos.size() << " GOs" << endl ;
	
	// write to file randset_out 
	std::ostream *out ;
	string outfile = directory + "/randset_out";
	out = new std::ofstream( outfile.c_str() ) ;
	
	*out << "Genes:\t" << gens.size() << endl ;
	*out << "GOs:\t" << gos.size() << endl ;
	
	// print line with names of GO-categories
	// print another line with sums of expressed genes annotated to the GOs (including background genes) 
	// also reset sums per GO to 0
	gos.print_names( *out ) ;

	
// 3) Read expressed test genes and add 1s to their annotated categories

	// read candidate_genes-file (infile data, one line per expressed TEST gene with genname
	string candidate_genes = directory + "/infile-data";
	std::ifstream candidate_genes_in( candidate_genes.c_str() ) ;
	if ( ! candidate_genes_in ) {
		Rcpp::stop("Cannot open infile-data.\n");
	}
	int n_candidate = 0 ; // number of expressed genes, previously size_of_random_sets
	string line_s ;	
	while( candidate_genes_in ) {
		getline( candidate_genes_in, line_s ) ;
		std::istringstream is( line_s.c_str() ) ;
		string gen_name ;
		is >> gen_name ;
		//Rcout << "expressed test gene: " << gen_name << endl;
		// genename_to_index: index of annotated GO-vec for a gene (gens[index])
		if ( genename_to_index.find( gen_name ) != genename_to_index.end() ) {
			// go through annotated GOs for the current gene and add 1
			for ( vector<int*>::iterator it = gens[genename_to_index[gen_name]].begin() ;
					it != gens[genename_to_index[gen_name]].end() ; ++it ) {
				(*(*it))++ ;
			}
			n_candidate++ ;
		}
	}
	candidate_genes_in.close();
	gos.print_sum( *out ) ;
	gos.clear() ;


// 4) Read regions and store in vector of bed-structures
	vector< bed_str > candidate_bed;
	vector< bed_str > background_bed;
	if (mod=="roll" || mod=="block"){	
		candidate_bed = read_bed(directory + "/test_regions.bed"); 
		background_bed = read_bed(directory+ "/bg_regions.bed"); 		
		// print regions 
		//Rcout << endl << "Candidate regions: " << endl;
		//for (int k=0; k < candidate_bed.size(); k++){
			//Rcout << candidate_bed[k].chrom << " " << candidate_bed[k].start << " " << candidate_bed[k].end  << " " << candidate_bed[k].len << " " << candidate_bed[k].cumu_len << endl;
		//}	
		//Rcout << endl << "Background regions: " << endl;
		//for (int k=0; k < background_bed.size(); k++){
			//Rcout << background_bed[k].chrom << " " << background_bed[k].start << " " << background_bed[k].end  << " " << background_bed[k].len << " " << background_bed[k].cumu_len << endl;
		//}	
	}


// 5) Create randomsets, randomly assign test genes among all expressed genes
	if (mod=="roll" || mod=="block"){
		Rcout << "Creating " << number_of_randomsets << " random gene sets from " << candidate_bed.size() << " random regions..." << endl ;	
	} else {
		Rcout << "Creating " << number_of_randomsets << " randomsets with size " << n_candidate << endl ;
	}			
	// for all randomsets
	for ( int i = 1 ; i <= number_of_randomsets ; ++i ) {
		set<int> random_numbers ; // indices for the selected genes
		
		if (mod=="block"){
			random_numbers = rannum_blocks(candidate_bed, background_bed, genename_to_index, genes_pos);
		} else if (mod=="roll"){
			random_numbers = rannum_roll(candidate_bed, background_bed, genename_to_index, genes_pos);
		} else if (mod=="gene_len"){			
			random_numbers = rannum_genelen(n_candidate, genename_to_index, genes_pos, total_length);
		} else if (mod=="classic"){
			while (random_numbers.size() < n_candidate) {	
				long ran = R::runif(0,1) * n_background;				
				random_numbers.insert(ran) ;
			}	
			//Rcout << "Candidate genes: " << n_candidate << ", Random genes: " << random_numbers.size() << std::endl;
		} 	
		
		// reset all go-nodes, set all go_vec counters to zero
		gos.clear(  ) ;

		// add 1 to every GO that a randomly choosen gene is part of
		// go through all random numbers 
		//Rcout << "random numbers: " << endl;
		for (set<int>::const_iterator it = random_numbers.begin() ;it != random_numbers.end() ; ++it){
			//Rcout << *it << ", ";
			// gens == genes. This vector is a simple representation of the go tree.
			// every gene is 1 vector of int*, where the int represents one go-node.
			// gens[*it] = vector with all GOs that belong to the gene *it points at
			for ( vector<int*>::iterator it2 = gens[*it].begin() ;it2 != gens[*it].end() ; ++it2 ) {
				(*(*it2))++ ;
			}
		}
		random_numbers.clear() ;
		
		// print a line with the sums of annotated genes for every go
		gos.print_sum( *out ) ;
	}
	delete out;
	genes_pos.clear();
	Rcout << "\rFinished" << endl ;
}
