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
