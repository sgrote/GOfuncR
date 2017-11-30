
#ifndef IDMAP_H
#define IDMAP_H

#include <cstdlib>
#include <iostream>
#include <string>
#include <map>
#define MAX_LINE_LENGTH_TERMS 200 

using std::string;
using std::istream;
using std::map;

/***************
 * Maps database ids (read from tab seperated file on construction)
 * to Go-Ids.
 *****************/
class idmap: public map<string, string>
{
	public:
		// term.txt from the 
		// go_termdb_tables distribution
		idmap( istream &in ) ;
		// get a database identifier for a GO-ID
		string get_id_for_go( string &go ) ;
	private:
} ;

#endif
