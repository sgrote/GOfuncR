
#include "go_graph_conti.h"
#include <sstream>

#include <Rcpp.h>
using std::endl;
using std::istringstream;

go_graph_conti::go_graph_conti( set<string> &nodes, istream &term2term, idmap &idm_ )
	: idm( idm_ ) 
{
	// read nodes, make go_obj_contis
	map<string, go_obj_conti*> temp_graph ;
	for ( set<string>::const_iterator it = nodes.begin() ; it!=nodes.end() ;
			++it ) 
	{
		temp_graph[*it] = new go_obj_conti( idm[*it] ) ;
	}

	// read term2term file, create graph
	while ( term2term ) {
		// 4 rows, 1 a id, 2 type of relationship, 3 parent, 4 child
		char s1[20] ;

		term2term.getline( s1, 20, '\t' ) ;
		term2term.getline( s1, 20, '\t' ) ;


		// parent id
		term2term.getline( s1, 20, '\t' ) ;

		map<string, go_obj_conti*>::const_iterator par = temp_graph.find( s1 ) ;
		if ( par != temp_graph.end() ) {
			// child id
			term2term.getline( s1, 20, '\n' ) ;
			string str_s1( s1 ) ;
			string child_id ;
			string::size_type tab_pos = str_s1.find( '\t' ) ;
			if ( tab_pos != string::npos ) {
				child_id = str_s1.substr( 0, tab_pos ) ;
			} else child_id = str_s1 ;
			map<string, go_obj_conti*>::const_iterator child = 
							temp_graph.find( child_id ) ;
			if ( child != temp_graph.end() ) {
				child->second->add_parent( par->second ) ;
			//	par->second->add_child( child->second ) ;
			}	
		}
	}
	// rewrite map file because detectedfile has
	// go_ids instead of database ids.
	for ( map<string, go_obj_conti*>::const_iterator it = temp_graph.begin() ;
			it != temp_graph.end() ; ++it ) 
	{
		graph[idm[it->first]] = it->second ;
	}
}

go_graph_conti::~go_graph_conti(  )
{
	for ( map<string, go_obj_conti*>::const_iterator it = graph.begin() ;
			it != graph.end() ; ++it ) 
	{
		delete it->second ;
	}
}

void go_graph_conti::get_parents( string &go, set<go_obj_conti*>* parents ) {
	
	if ( graph.find( go ) != graph.end() ) {
		graph[go]->get_parents( parents ) ;
	} else {
		Rcpp::Rcerr << "Cannot find " << go << endl ;
	} 
}


void go_graph_conti::clear_vals(  ) {
	for ( map<string, go_obj_conti*>::const_iterator it = graph.begin() ;
			it != graph.end() ; ++it ) 
	{
		it->second->clear( ) ;
	}
}

void go_graph_conti::print_groups( ostream &os ) 
{
	for ( map<string,go_obj_conti*>::const_iterator it = graph.begin() ;
				it != graph.end() ; ++it ) 
	{
		os << it->first << '\t' ;
	}
	os << '\n' ;
}

void go_graph_conti::print_vals( ostream &os )
{
	for ( map<string,go_obj_conti*>::const_iterator it = graph.begin() ;
				it != graph.end() ; ++it ) 
	{
		it->second->print( os ) ; 
	}
	os << '\n' ;
	
}

void go_graph_conti::print_nr_genes( ostream &os )
{
	for ( map<string,go_obj_conti*>::const_iterator it = graph.begin() ;
				it != graph.end() ; ++it ) 
	{
		it->second->print_nr_genes( os ) ; 
	}
}
