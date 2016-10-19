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


struct gen_pos_str {
	string name; string chrom; int start; int end; 
};	
struct bed_str {
	string chrom; int start; int end; double len_quant;
};	

//[[Rcpp::export]]
void hyper_randset_blocks(std::string detected, std::string changed ,int number_of_randomsets, std::string outfile, std::string term, std::string graph_path, std::string termtoterm, std::string root) 

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
	


// 2) Read expressed genes and their annotated GO-categories
		// gos-object will be used to get one int* per GO and to print the results.
	go gos ;	
	std::ostream *out ;
	out = new std::ofstream( outfile.c_str() ) ;
	std::ifstream in( detected.c_str() ) ;
	Rcerr << "Reading detectedfile... " << endl ;
			
	// gens == genes. This vector is a simple representation of the go tree.
	// every gene is 1 vector of int*, where the int* represents one go-node.
	vector< vector<int*> > gens;
	// vector of structures to store coordinates of genes from Allen:4005
	vector<gen_pos_str> genes_pos;
	// genename to index for gens
	map<string,int> genename_to_index ;
	int index = 0 ;

	string line ;
	// in = detected = Allen:4005-file = one line per expr gene with genename and all its anno. GOs
	while ( getline(in, line) ) { 
		std::istringstream is( line.c_str() ) ; // take line as input stream
		string gen_name ;
		is >> gen_name ; // speichere erstes wort in line (genename) in gen_name			
#ifdef DEBUG		
		Rcout << "gen_name: " << gen_name << endl ;
#endif
		// store position of current gene
		gen_pos_str gen_pos; 
		gen_pos.name = gen_name;
		is >> gen_pos.chrom;
		is >> gen_pos.start;
		is >> gen_pos.end;
		genes_pos.push_back(gen_pos);
		
		// store pointers to GOs for the current gene
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
		Rcout << "expressed test gene: " << gen_name << endl;
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

	
// 3) Read desert regions and get length of deserts 

	vector<int> desert_lens;
	int chrom, start, end, length;	
	std::ifstream desert ("/r1/people/steffi_grote/ownCloud/forAkeyPaper/deserts.bed");
	while(getline(desert, line)){
		std::istringstream is( line.c_str() ); 
		is >> chrom; is >> start; is >> end;
		length = end - start;
		desert_lens.push_back(length);
	}
	Rcout << "Desert region lengths:" << endl;
	for (int i=0; i<desert_lens.size(); i++){
		Rcout << desert_lens[i] << endl; 
	}
	
	
//// 4) Read mappable regions and store in vector of bed-structures

	vector< bed_str > mappable_bed;
	std::ifstream mappable ("/r1/people/steffi_grote/ownCloud/forAkeyPaper/mappable_regions.bed");
	while(getline( mappable, line )){
		bed_str bed;
		std::istringstream is( line.c_str() ) ;
		is >> bed.chrom; is >> bed.start; is >> bed.end;		
		mappable_bed.push_back(bed);  
	}


// 5) Create randomsets, randomly assign test genes among all expressed genes
	Rcerr << "Creating " << number_of_randomsets << " random gene sets from " << desert_lens.size() << " random deserts..." << endl ;	
			
	// forall randomsets
	for ( int i = 1 ; i <= number_of_randomsets ; ++i ) {
		
		// create randomset
		set<string> random_genes; // names of selected genes
		set<int> random_numbers ; // indices for the selected genes
		int sum_genes = 0;

		// loop over deserts
		for (int j=0; j < desert_lens.size(); j++){
		// for (int j=0; j < 2; j++){
			// randomly choose a mappable region
			// random [0,1] to choose from quantiles of total mappable length
			int tot_mappable_len = 0; 
			for (int k=0; k < mappable_bed.size(); k++){  // get total length of mappable regions 
				int this_len = mappable_bed[k].end - mappable_bed[k].start - desert_lens[j];
				if (this_len < 0) this_len = 0;
				tot_mappable_len += this_len;
			}
			
			int cumu_len = 0;
			for (int k=0; k < mappable_bed.size(); k++){  // get length quantiles for each region
				int this_len = mappable_bed[k].end - mappable_bed[k].start - desert_lens[j];
				if (this_len < 0) this_len = 0;
				cumu_len += this_len;
				mappable_bed[k].len_quant = cumu_len*1.0 / tot_mappable_len ;
			}
			
			//double ran = (double) rand() / (RAND_MAX);
			double ran = R::runif(0,1);
			int k = 0;
			while (mappable_bed[k].len_quant < ran){  // choose chromosome and region
				k++;
			}	
			int begin_start = mappable_bed[k].start;
			int end_start = mappable_bed[k].end - desert_lens[j] ;
			string ran_chrom = mappable_bed[k].chrom; 
			
			// randomly choose a starting point between begin_start and end_start
			//double ran2 = (double) rand() / (RAND_MAX);
			double ran2 = R::runif(0,1);
			int ran_start = ran2 * (end_start - begin_start) + begin_start;
			int ran_end = ran_start + desert_lens[j];
			
			//Rcout << endl << ran_chrom << " " << ran_start << " " << ran_end << " " << ran_end-ran_start << endl;
		
	
// 6) go through genes positions and select those that overlap randomly chosen region
			// genes positions stored in genes_pos<genes_pos_str>
			for (int g=0; g<genes_pos.size(); g++){
				if(genes_pos[g].chrom == ran_chrom &&
				((genes_pos[g].start > ran_start && genes_pos[g].start < ran_end) ||
				(genes_pos[g].end > ran_start && genes_pos[g].end < ran_end))){
					// add to set of randomly chosen test genes	
					random_genes.insert(genes_pos[g].name);  
					sum_genes ++;
					//Rcout << genes_pos[g].name << " " << genes_pos[g].chrom << " " << genes_pos[g].start << " " << genes_pos[g].end << " " << genename_to_index[genes_pos[g].name] << endl;			
				}
			}
		}
		Rcout << "sum of random genes: " << sum_genes << endl;
		Rcout << "sum of unique random genes: " << random_genes.size() << endl;

	//cout << "Mappable regions:" << endl;
	//for (int i=0; i<mappable_bed.size(); i++){
		//cout << mappable_bed[i].chrom << " " << mappable_bed[i].start << " " << mappable_bed[i].end << " " << mappable_bed[i].len_quant << endl;
	//}
	
// 7) sum counts for ontology terms annotated to the randomly selected genes   
		// get index for every random gene name
		for (set<string>::iterator it=random_genes.begin(); it!=random_genes.end(); ++it){
			random_numbers.insert(genename_to_index[*it]);
		}
		//Rcout << endl << "random_genes and their index for 'gens' with GO-vec:" << endl;
		//for (set<string>::iterator it=random_genes.begin(); it!=random_genes.end(); ++it){
			//Rcout << *it << " " << genename_to_index[*it] << endl;
		//}
		
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
		random_genes.clear();
		random_numbers.clear() ;
		
		// print a line with the sums of annotated genes for every go
		gos.print_sum( *out ) ;
	}
	delete out;
	desert_lens.clear();
	mappable_bed.clear();
	genes_pos.clear();
	Rcerr << "\rFinished" << endl ;
}
