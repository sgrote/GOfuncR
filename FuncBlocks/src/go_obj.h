
#ifndef GO_OBJ_H
#define GO_OBJ_H

#include <string>
#include <vector>
#include <set>
#include <iostream>

using std::vector ;
using std::set ;
using std::ostream ;
using std::string ;

class gene ;


/************
 * capsulates information a GO-Group contains:
 * Name, Number Genes, rank for wilcoxon test
 *************/
class go_obj
{
	public:
		go_obj( string &name_ ) ;
		void add_parent( go_obj* p ) ;
		void get_parents( set<go_obj*> *parentset ) ;
		void add_gene( gene* g_ ) ; 

		// genes are swapping annotation to 
		// create random distribution, clear before new
		// genes annotate itself
		void clear_genes(  ) ; 

		// sum of ranks of all annotated genes
		void print_sumranks( ostream &os ) ;

		// number genes registered in this node
		void print_n( ostream &os ) ;

	private:
		string name ;
		vector<go_obj*> parents ;
		vector<gene*> genes ;
		
} ;

#endif
