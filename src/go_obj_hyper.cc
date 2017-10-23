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


#include "go_obj_hyper.h" 
#include <cmath>

using std::endl;

go_obj_hyper::go_obj_hyper( string &name_ ) : name( name_ )
{  }

void go_obj_hyper::add_parent( go_obj_hyper* p )
{ parents.push_back( p ) ; }

void go_obj_hyper::get_parents( set<string> *parentset ) {
	parentset->insert( parentset->begin(), name ) ;
	for ( vector<go_obj_hyper*>::const_iterator it = parents.begin() ;
			it != parents.end() ; ++it ) 
	{
		(*it)->get_parents( parentset ) ;
	}
}

