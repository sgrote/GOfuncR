
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

