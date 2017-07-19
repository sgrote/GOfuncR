
#ifndef GO_GRAPH_CONTI_H
#define GO_GRAPH_CONTI_H

#include "go_obj_conti.h"
#include "idmap.h"
#include <map>
#include <set>
#include <iostream>

using std::ostream;

/******
 * Reads GO DAG and handles go_objs. 
 *******/
class go_graph_conti
{
	public:
		/*********
 		 * Reads DAG of GO and creates go_objs. 
		 * idmap: map database id of termdb-tables to GO_Id
		 * nodes: all nodes below the root node
		 * term2term: file with term2term associations
		 ***********/
		go_graph_conti( set<string> &nodes, istream &term2term, idmap &idm_ ) ;
		~go_graph_conti(  ) ;

		/*********
 		 * returns parent go_objs of node with name "go", 
		 * including the node itself
		 ***********/
		void get_parents( string &go, set<go_obj_conti*> *parents ) ;

		/*********
 		 * delete data of all go_objs
		 ***********/
		void clear_vals(  ) ;

		/*********
 		 * name of groups, values, separated by tab 
		 * on one line
		 * (graph is not touched, so order remains the same 
		 * for all 2 functions)
		 ***********/
		void print_groups( ostream &os ) ;
		void print_vals( ostream &os ) ;

		// number of genes for each go. Format: GO \t #genes\n
		void print_nr_genes( ostream &os ) ;

	private:
		idmap &idm ;
		map<string, go_obj_conti*> graph ;
} ;


#endif
