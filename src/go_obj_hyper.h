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




#ifndef GO_OBJ_HYPE_H
#define GO_OBJ_HYPE_H

#include <string>
#include <vector>
#include <set>
#include <iostream>

using std::vector ;
using std::set ;
using std::ostream ;
using std::string ;

/**
 Represents a single GO node.
 */
class go_obj_hyper
{
	public:
		// name_ is a identifier (i.e. GO:Id). 
		go_obj_hyper( string &name_ ) ;
		// p is a direct parent of this GO-Node
		void add_parent( go_obj_hyper* p ) ;
		// get all parents that are direct or indirect parents of this GO-Node.
		void get_parents( set<string> *parentset ) ;

	private:
		string name ;
		vector<go_obj_hyper*> parents ;
} ;

#endif
