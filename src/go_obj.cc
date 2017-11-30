
#include "go_obj.h" 
#include <math.h>
#include "gene.h"

#include <Rcpp.h>

using std::endl;

go_obj::go_obj( string &name_ ) : name( name_ )
{  }

void go_obj::add_parent( go_obj* p )
{ parents.push_back( p ) ; }

void go_obj::get_parents( set<go_obj*> *parentset ) {
	parentset->insert( parentset->begin(), this ) ;
	for ( vector<go_obj*>::const_iterator it = parents.begin() ;
			it != parents.end() ; ++it ) 
	{
		(*it)->get_parents( parentset ) ;
	}
}


void go_obj::add_gene( gene* g_ ) 
{ 
	genes.push_back( g_ ) ;
} 

void go_obj::clear_genes(  ) 
{
	genes.clear(  ) ;
}

void go_obj::print_sumranks( ostream &os ) 
{ 
	double sumranks = 0. ;
	for ( vector<gene*>::const_iterator it = genes.begin() ;
			it != genes.end() ; ++it ) 
	{
		if ( (*it)->get_rank() == 0 ) {
			Rcpp::Rcerr << name << ": " << (*it)->name << " rank == 0 " << endl ;
		}
		sumranks += (*it)->get_rank() ;
	}  
	os << sumranks << '\t' ;
}

void go_obj::print_n( ostream &os ) 
{ 
	os << genes.size() << '\t' ; 
}

