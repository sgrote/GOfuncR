
#ifndef GENES_CONTI_H
#define GENES_CONTI_H

#include <fstream>
#include <map>
#include <sstream>
#include "go_graph_conti.h"
#include "gene_conti.h"
#include <list>

using namespace std ;

/********************
 * handles gene-objects and shuffle gene_data
 *********************/
class genes_conti
{
	public:
		// annotationfile = Gene \t GO_1\ GO_2 ...
		// datafile = Gene \t leftvalue \t rightvalue
		genes_conti( go_graph_conti &graph, istream &annotation, istream &data ) ;

		// shuffle genes and exchange data of genes.
		void create_random_set(  ) ;
		~genes_conti(  ) ;
	private:
		map<string, gene_conti*> genemap ;
		vector<gene_conti*> genevec ;

} ;

#endif
