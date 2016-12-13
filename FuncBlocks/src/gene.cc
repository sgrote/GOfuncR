
#include "gene.h"
#include <time.h>
#include <cstdlib>
#include <cstdio>

void gene::write_to_gos( set<go_obj*>* gos_ ) 
{
	for ( set<go_obj*>::const_iterator it = gos_->begin() ;
			it != gos_->end() ; ++it ) 
	{
		(*it)->add_gene( this ) ;
	}
	
}

void gene::write_to_gos(  ) 
{
	for ( set<go_obj*>::const_iterator it = gos.begin() ; it != gos.end() ; ++it ) {
		(*it)->add_gene( this ) ;
	}
	
}

void gene::set_rank( double _rank ) 
{
	rank = _rank ;
}

double gene::get_rank(  ) 
{
	return( rank ) ;
}

set<go_obj*>* gene::get_gos( ) 
{
	return &gos ;
}

