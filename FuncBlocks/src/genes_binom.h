
#ifndef GENES_BINOM_H
#define GENES_BINOM_H

#include <fstream>
#include <map>
#include <sstream>
#include "go_graph_binom.h"
#include "gene_binom.h"

/********************
 * handles gene-objects and shuffle gene_data
 *********************/
class genes_binom
{
	public:
		// annotationfile = Gene \t GO_1\ GO_2 ...
		// datafile = Gene \t leftvalue \t rightvalue
		genes_binom( go_graph_binom &graph, istream &annotation, istream &data ) ;

		// shuffle genes and exchange data of genes.
		void create_random_set(  ) ;
		~genes_binom(  ) ;
	private:
		map<string, gene_binom*> genemap ;
		vector<gene_binom*> gene_vec ;

} ;

#endif
