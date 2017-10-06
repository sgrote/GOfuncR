
// verhindert dass go_graph merhmals included wird
#ifndef GO_GRAPH_H
#define GO_GRAPH_H

#include "go_obj.h"
#include "idmap.h"
#include <map>
#include <set>
#include <iostream>
#include <string>

using std::ostream;

/******
 * Reads GO DAG and handles go_objs. 
 *******/
class go_graph
{
	public:
		/*********
 		 * Reads DAG of GO and creates go_objs. 
		 * idmap: map database id of termdb-tables to GO_Id
		 * nodes: all nodes below the root node
		 * term2term: file with term2term associations
		 ***********/
		go_graph( set<string> &nodes, istream &term2term, idmap &idm_ ) ;
		~go_graph(  ) ;

		/*********
 		 * returns parent go_objs of node with name "go", 
		 * including the node itself
		 ***********/
		void get_parents( string &go, set<go_obj*> *parents ) ;

		/*********
 		 * instead of swapping data between genes, simply annotation 
		 * is swapped. clear_genes kills annotation in all nodes. see
		 * genes.h: shuffle_genes
		 ***********/
		void clear_genes(  ) ;

		
		/*********
 		 * header(): name of groups, number genes separated by tab 
		 * sumranks(): sum of all ranks in groups seperated by tab
		 * (graph is not touched, so order remains the same 
		 * for all 2 functions)
		 ***********/
		void print_header( ostream &os ) ;
		void print_sumranks( ostream &os ) ;
	private:
		idmap &idm ;
		map<string, go_obj*> graph ;
} ;


#endif
