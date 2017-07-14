
#include "go_obj_binom.h" 
#include <math.h>

using std::cerr;
using std::cout;
using std::endl;

go_obj_binom::go_obj_binom( string &name_ ) : name( name_ ), hka(0), cka(0), nr_genes(0)
{  }

void go_obj_binom::add_parent( go_obj_binom* p )
{ parents.push_back( p ) ; }

void go_obj_binom::get_parents( set<go_obj_binom*> *parentset ) {
	parentset->insert( parentset->begin(), this ) ;
	for ( vector<go_obj_binom*>::const_iterator it = parents.begin() ;
			it != parents.end() ; ++it ) 
	{
		(*it)->get_parents( parentset ) ;
	}
}


void go_obj_binom::add_hka( int x=1 ) 
{ 
	hka+=x ; 
} 

void go_obj_binom::add_cka( int x=1 ) 
{ 
	cka+=x ; 
}

void go_obj_binom::add_gene() 
{
	nr_genes++ ;
}

void go_obj_binom::clear_ka(  ) 
{ 
	hka = 0 ; 
	cka = 0 ; 
} 

void go_obj_binom::print_nr_genes( ostream &os ) 
{
	os << name << '\t' << nr_genes << endl ;
}

void go_obj_binom::print_ka( ostream &os ) 
{ 
	os << hka << ' ' << cka << '\t' ; 
}

