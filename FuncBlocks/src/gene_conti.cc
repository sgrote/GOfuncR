
#include "gene_conti.h"

void gene_conti::add( int ch_s_, int ch_ns_, int hh_s_, int hh_ns_ ) {
	ch_s = ch_s_ ;
	ch_ns = ch_ns_ ;
	hh_s = hh_s_ ;
	hh_ns = hh_ns_ ;
	for ( set<go_obj_conti*>::const_iterator it = gos.begin() ;
			it != gos.end() ; ++it ) 
	{
		(*it)->add( ch_s, ch_ns, hh_s, hh_ns ) ;
		(*it)->add_gene(  ) ;
	}
}

void gene_conti::write_data_gos( set<go_obj_conti*>* gos_ ) 
{
	for ( set<go_obj_conti*>::const_iterator it = gos_->begin() ;
			it != gos_->end() ; ++it ) 
	{
		(*it)->add( ch_s, ch_ns, hh_s, hh_ns ) ;
	}
	
}

set<go_obj_conti*>* gene_conti::get_gos( ) 
{
	return &gos ;
}

