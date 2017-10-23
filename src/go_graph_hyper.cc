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




#include "go_graph_hyper.h"
#include <sstream>

#include <Rcpp.h>

using std::endl;
using std::istringstream;

go_graph_hyper::go_graph_hyper( set<string> &nodes, istream &term2term, idmap &idm_ )
	: idm( idm_ ) 
{
	// read nodes, make go_objs
	map<string, go_obj_hyper*> temp_graph ;
	for ( set<string>::const_iterator it = nodes.begin() ; it!=nodes.end() ;
			++it ) 
	{
		temp_graph[*it] = new go_obj_hyper( idm[*it] ) ;
	}

	// read term2term file, create graph
	while ( term2term ) {
		// 4 rows, 1 a id, 2 type of relationship, 3 parent, 4 child
		char s1[20] ;

		term2term.getline( s1, 20, '\t' ) ;
		term2term.getline( s1, 20, '\t' ) ;


		// parent id
		term2term.getline( s1, 20, '\t' ) ;

		map<string, go_obj_hyper*>::const_iterator par = temp_graph.find( s1 ) ;
		if ( par != temp_graph.end() ) {
			// child id
			term2term.getline( s1, 20, '\n' ) ;
			string str_s1( s1 ) ;
			string child_id ;
			string::size_type tab_pos = str_s1.find( '\t' ) ;
			if ( tab_pos != string::npos ) {
				child_id = str_s1.substr( 0, tab_pos ) ;
			} else child_id = str_s1 ;
			map<string, go_obj_hyper*>::const_iterator child = 
												temp_graph.find( child_id ) ;
			if ( child != temp_graph.end() ) {
				child->second->add_parent( par->second ) ;
			//	par->second->add_child( child->second ) ;
			}	
		} else {
			// skip rest if the parent node is not part of the graph
			term2term.getline( s1, 20, '\n' ) ;
		}
	}
	// rewrite map file because detectedfile has
	// go_ids instead of database ids.
	for ( map<string, go_obj_hyper*>::const_iterator it = temp_graph.begin() ;
			it != temp_graph.end() ; ++it ) 
	{
		graph[idm[it->first]] = it->second ;
	}
}

go_graph_hyper::~go_graph_hyper(  )
{
	for ( map<string, go_obj_hyper*>::const_iterator it = graph.begin() ;
			it != graph.end() ; ++it ) 
	{
		delete it->second ;
	}
}

void go_graph_hyper::get_parents( string &go, set<string>* parents ) {
	
	if ( graph[go] ) {
		graph[go]->get_parents( parents ) ;
	} else {
		Rcpp::Rcerr << "Error: Cannot find " << go 
			<< ". Maybe taxonomies are not in the right order." << endl ;
	}
}
