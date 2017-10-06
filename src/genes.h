
#ifndef GENES_H
#define GENES_H

#include <fstream>
#include <map>
#include <sstream>
#include "go_graph.h"
#include "gene.h"

/*********
 * handles gene objects
 **********/
class genes
{
	public:
		
		/*********
		 * annotationfile = Gene \t GO_1\ GO_2 ...
		 * datafile = Gene \t float
		 * read data, create gene-objects with data and annotation
		 **********/
		genes( go_graph &graph, istream &annotation, istream &data ) ;

		// needed for calculation of wilcoxon rank
		double sumnties( ) { return sum_nties ; } ;

		// shuffle gene annotation
		void create_random_set(  ) ;
		~genes(  ) ;
	private:
		double sum_nties ;
		map<string, gene*> genemap ;
		vector<gene*> gene_vec ;

} ;

#endif
