
#include "go_obj.h"

using namespace std ;

/*****************
 * gene encapsulates relevant gene information: name, GO annotation, data
 ******************/
class gene {
	public:
		
		gene( string name_, set<go_obj*> &gos_ ) : name( name_ ), gos(gos_)
		{  }
		set<go_obj*>* get_gos(  ) ;

		// write rank to gos_
		void write_to_gos( set<go_obj*>* gos_ ) ;

		// write rank to private: gos 
		void write_to_gos(  ) ;

		void set_rank( double rank_ ) ;
		double get_rank(  ) ;
		
		string name ;
	private:
		set<go_obj*> gos ;
		double rank ;

} ;
