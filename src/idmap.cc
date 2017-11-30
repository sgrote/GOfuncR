
#include "idmap.h"
#include <sstream>

idmap::idmap( istream &in ) 
{
	char line[MAX_LINE_LENGTH_TERMS] ;

	while ( in ) {
		// first row: id
		in.getline( line, MAX_LINE_LENGTH_TERMS, '\t' ) ;
		string id( line ) ;
		
		// skip 2 fields
		in.getline( line, MAX_LINE_LENGTH_TERMS, '\t' ) ;
		in.getline( line, MAX_LINE_LENGTH_TERMS, '\t' ) ;

		// GO:Number
		in.getline( line, MAX_LINE_LENGTH_TERMS, '\t' ) ;
		string go( line ) ;
		
		if ( id.size() > 0 && go.size() > 0 ) (*this)[id]=go ;
		
		in.getline( line, MAX_LINE_LENGTH_TERMS, '\n' ) ;
	}
}

string idmap::get_id_for_go( string &go ) 
{
	for ( map<string,string>::const_iterator it = this->begin() ; 
			it != this->end() ; ++it ) 
			if ( it->second == go ) return it->first ;

	return 0 ;
}
