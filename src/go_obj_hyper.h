
#ifndef GO_OBJ_HYPE_H
#define GO_OBJ_HYPE_H

#include <string>
#include <vector>
#include <set>
#include <iostream>

using std::vector ;
using std::set ;
using std::ostream ;
using std::string ;

/**
 Represents a single GO node.
 */
class go_obj_hyper
{
	public:
		// name_ is a identifier (i.e. GO:Id). 
		go_obj_hyper( string &name_ ) ;
		// p is a direct parent of this GO-Node
		void add_parent( go_obj_hyper* p ) ;
		// get all parents that are direct or indirect parents of this GO-Node.
		void get_parents( set<string> *parentset ) ;

	private:
		string name ;
		vector<go_obj_hyper*> parents ;
} ;

#endif
