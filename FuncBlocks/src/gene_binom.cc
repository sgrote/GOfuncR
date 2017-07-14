
#include "gene_binom.h"
//#include <time.h>
#include <cstdlib>
#include <cstdio>

void gene_binom::add_hka_cka( int hka_, int cka_ ) {
	
	hka = hka_ ;
	cka = cka_ ;

	for ( set<go_obj_binom*>::const_iterator it = gos.begin() ;
			it != gos.end() ; ++it ) 
	{
		(*it)->add_cka( cka ) ;
		(*it)->add_hka( hka ) ;
		(*it)->add_gene(  ) ;
	}
}

void gene_binom::write_data_gos( set<go_obj_binom*>* gos_ ) 
{
	for ( set<go_obj_binom*>::const_iterator it = gos_->begin() ;
			it != gos_->end() ; ++it ) 
	{
		(*it)->add_cka( cka ) ;
		(*it)->add_hka( hka ) ;
	}
	
}

set<go_obj_binom*>* gene_binom::get_gos( ) 
{
	return &gos ;
}

