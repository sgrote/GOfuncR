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
		// in should be the file terms.txt from the 
		// go_termdb_tables distribution
		idmap( istream &in ) ;
		// get a database identifier for a GO-ID
		string get_id_for_go( string &go ) ;
	private:
} ;

#endif
