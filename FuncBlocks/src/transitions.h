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

