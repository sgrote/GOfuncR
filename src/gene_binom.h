
#ifndef GENE_BINOM_H
#define GENE_BINOM_H


#include "go_obj_binom.h"

using namespace std ;

/***********
 * encapsulates information of a single gene
 ************/
class gene_binom {
	public:
		/********
                 * name_ == name of gene
		 * gos == groups gene is annotated to (including subsumtion)
		 *********/
		gene_binom( string name_, set<go_obj_binom*> &gos_ ) : name( name_ ), gos(gos_)
		{  }


		/**********
                 * adds data to GO groups, including information 
		 * that a gene is annotated to the group and save data to hka,cka
		 ***********/
		void add_hka_cka( int hka_, int cka_ ) ;

		set<go_obj_binom*>* get_gos(  ) ; 

		/**********
		 * adds data to GO groups
		 ***********/
		void write_data_gos( set<go_obj_binom*>* gos_ ) ;
		
		string name ;
	private:
		set<go_obj_binom*> gos ;
		int hka, cka ; // right and left value for binomial test

} ;

#endif
