
#ifndef TRANSI_H
#define TRANSI_H

#include <string>
#include <set>
#include <iostream>

using std::string;
using std::set;
using std::istream;

#define MAX_LINE_LENGTH_TRANS 100 

/***********************
 * Given a GO-ID and the graph_path.txt file from the termdb_tables distribution of
 * Gene Ontology an object with all GO-Ids will be created, that are reachable from
 * "id".
 ******************************/
class transitions: public set<string>
{
	public:
		transitions( string &id, istream &in ) ;
} ;

#endif
