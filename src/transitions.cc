
#include "transitions.h"

transitions::transitions( string &id, istream &in ) 
{
	char line[MAX_LINE_LENGTH_TRANS] ;

	while( in ) {
		// first row, id, skip
		in.getline( line, MAX_LINE_LENGTH_TRANS, '\t' ) ;

		// second line, parent
		in.getline( line, MAX_LINE_LENGTH_TRANS, '\t' ) ;
		string parent( line ) ;
		if ( parent == id ) {
			in.getline( line, MAX_LINE_LENGTH_TRANS, '\t' ) ;
			string child( line ) ;
			this->insert( this->begin(), child ) ;
		}
		// skip
		in.getline( line, MAX_LINE_LENGTH_TRANS, '\n' ) ;
	}


	// also add parent
	this->insert( id ) ;

}
