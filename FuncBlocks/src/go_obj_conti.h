
#ifndef GO_OBJ_CONTI_H
#define GO_OBJ_CONTI_H

#include <string>
#include <vector>
#include <set>
#include <iostream>

using std::vector ;
using std::set ;
using std::ostream ;
using std::string ;

/************
 * encapsulates information a GO-Group contains:
 * Name, Number Genes, left and right values for binom test 
 *************/
class go_obj_conti
{
	public:
		go_obj_conti( string &name_ ) ;
		void add_parent( go_obj_conti* p ) ;
		void get_parents( set<go_obj_conti*> *parentset ) ;
		void add( int ch_s_, int ch_ns_, int hh_s_, int hh_ns_ ) ;
		void add_gene( ) ; 
		void clear(  ) ; 
		void print( ostream &os ) ;
		void print_nr_genes( ostream &os ) ;

	private:
		string name ;
		vector<go_obj_conti*> parents ;
		int ch_s, ch_ns, hh_s, hh_ns ;
		int nr_genes ;
} ;

#endif
