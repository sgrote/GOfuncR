/*
 * FUNC - Functional Analysis of Gene Expression Data
 * Copyright (C) 2002  Bjoern Muetzel, Kay Pruefer
 * 
 * This program is modifiable/redistributable under the terms of the
 * GNU General Public License.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */




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
