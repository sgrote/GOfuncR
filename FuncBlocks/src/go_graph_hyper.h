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




#ifndef GO_GRAPH_HYPE_H
#define GO_GRAPH_HYPE_H

#include "go_obj_hyper.h"
//steffi:
#include "idmap.h"
#include <map>
#include <set>

using std::ostream;

/*
 * a class to handle go_objs
 */
class go_graph_hyper
{
	public:
		// "noded" are the nodes which will be used to build the graph
		// term2term == term2term.txt from termdb_tables distribution of GO
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
