
#include "go.h"

int *go::add( const string &go_name ) 
{
	map<string, int*>::iterator it = go_map.find( go_name ) ;
	if ( it != go_map.end() ) { // found
		(*(it->second))++ ;
		return it->second ;
	} else {
		int *x = new int ;
		*x = 1 ;
		go_map[go_name] = x ;
		return x ;
	}
}

void go::clear(  ) 
{
	for ( vector<int*>::iterator it = go_vec.begin() ; it != go_vec.end() ;
			++it ) 
		**it = 0 ;
}

void go::print_sum( ostream &os ) 
{
	for ( vector<int*>::iterator it = go_vec.begin() ; it != go_vec.end() ;
			++it ) 
		os << **it << "\t" ;
	os << endl ;
}
void go::print_names( ostream &os ) 
{
	for ( map<string,int*>::iterator it = go_map.begin() ; it != go_map.end() ; )
	{
		os << it->first << "\t" ;
		go_vec.push_back( it->second ) ;
		map<string,int*>::iterator it2 = it ;
		++it2 ;
		go_map.erase( it ) ;
		it = it2 ;
	}
	os << endl ;
	for ( vector<int*>::iterator it = go_vec.begin() ; it != go_vec.end() ;
			++it ) 
	{
		os << **it << "\t" ;
		**it = 0 ;
	}
	os << endl ;
}

go::~go(  ) 
{
	for ( vector<int*>::const_iterator it = go_vec.begin() ; it != go_vec.end() ; ++it )
		delete *it ;
	
}

