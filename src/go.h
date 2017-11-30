
#ifndef GO_H
#define GO_H


#include <map>
#include <vector>
#include <string>
#include <iostream>

using std::map;
using std::vector;
using std::string;
using std::ostream;
using std::endl;


/*
 * Class is used to print GO values (howmany changed genes per GO),
 * reset these values and to get a unique int* for every GO.
 */
class go {
	public:
		go() {} 

		// add a new go_name and return pointer to counter
		int *add( const string &go_name ) ;
		
		// set all counters to zero
	    void clear(  ) ;

		// print sum for each counter (tab seperated)
		void print_sum( ostream &os ) ;

		// print all go-names (or ids) (tab seperated)
		// print the number of genes for every go id
		// side_effekt: go_map will be converted to go_vec
		void print_names( ostream &os ) ;

		// nr of GO-Terms.
		int size() { return go_map.size() ; } ;
		
		// kill all int* allocated somewhere else
		~go(  ) ;
	private:
		map<string,int*> go_map ;
		vector<int*> go_vec ;
} ;

#endif
