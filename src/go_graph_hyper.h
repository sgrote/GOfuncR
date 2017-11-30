
#ifndef GO_GRAPH_HYPE_H
#define GO_GRAPH_HYPE_H

#include "go_obj_hyper.h"
#include "idmap.h"
#include <map>
#include <set>

using std::ostream;

// a class to handle go_objs

class go_graph_hyper
{
	public:
		// idmap: see idmap.h
		go_graph_hyper( set<string> &nodes, istream &term2term, idmap &idm_ ) ;
		~go_graph_hyper(  ) ;
		// get all (direct or indirect) parents for GO-ID go
		void get_parents( string &go, set<string> *parents ) ;
	private:
		idmap &idm ;
		map<string, go_obj_hyper*> graph ;
} ;


#endif
