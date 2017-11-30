
#ifndef FISH_TEST_H
#define FISH_TEST_H

#include <cmath>
#include <iostream> 
#include <Rcpp.h>
using namespace std ;

double logfak( int x ) {
	if ( x <= 1 ) return 0. ;
	return log(static_cast<double>(x))+logfak(x-1) ;
}

///////////////////////////////////////////////////////////////
// FISHERS EXACT TEST, 2-tail
///////////////////////////////////////////////////////////////
//
// Calculates fishers exact with an 2x2 table, 2-tail
// Arguments:  |   |   |
// 	       +--------
//             | a1| b1|
//             | a2| b2|
// See: Sokahl, Rohlf "Biometry" p.734
///////////////////////////////////////////////////////////////
	
double disproportion( int a1, int a2, int b1, int b2 ) 
{
	return( abs( static_cast<double>(a1)/static_cast<double>(a1+b1)
					 - static_cast<double>(a2)/static_cast<double>(a2+b2) ) ) ;

}


double fishers_exact_2t( int a1, int a2, int b1, int b2 ) {
	double sum = logfak( a1+a2 ) + logfak( a1+b1 ) + logfak( b1+b2 ) + logfak( a2+b2 ) - logfak( a1+b1+a2+b2 ) ;
	
	double data_disprop = disproportion( a1, a2, b1, b2 ) ;
	double prob = 0.0 ;
	int orig_a1 = a1 ;

	if ( (a1*b2 - a2*b1) < 0 ) {
		while ( a1 >=0 && b2 >= 0 ) {
			prob += exp( sum - ( logfak(a1)+logfak(a2)+logfak(b1)+logfak(b2) ) ) ;

			a1-- ;
			b2-- ;
			a2++ ;
			b1++ ;
		}

		a1++; b2++; a2--; b1--;

		if ( a2 < b1 ) {
			a1+=a2 ;
			b2+=a2 ;
			b1-=a2 ;
			a2 = 0 ;
		} else {
			a1+=b1 ;
			b2+=b1 ;
			a2-=b1 ;
			b1 = 0 ;
		}
		while( disproportion( a1, a2, b1, b2 ) > data_disprop ) {
			prob += exp( sum - ( logfak(a1)+logfak(a2)+logfak(b1)+logfak(b2) ) ) ;
			a1-- ;
			b2-- ;
			a2++ ;
			b1++ ;
			if ( a1 == orig_a1 ) break ;
		}
	} else {
		while ( a2 >=0 && b1 >= 0 ) {
			prob += exp( sum - ( logfak(a1)+logfak(a2)+logfak(b1)+logfak(b2) ) ) ;
			a1++ ;
			b2++ ;
			a2-- ;
			b1-- ;
		}
		a1--; b2--; a2++; b1++;

		if ( a1 < b2 ) {
			a2+=a1 ;
			b1+=a1 ;
			b2-=a1 ;
			a1 = 0 ;
		} else {
			a2+=b2 ;
			b1+=b2 ;
			a1-=b2 ;
			b2 = 0 ;
		}

		while( disproportion( a1, a2, b1, b2 ) > data_disprop ) {
			prob += exp( sum - ( logfak(a1)+logfak(a2)+logfak(b1)+logfak(b2) ) ) ;
			a1++ ;
			b2++ ;
			a2-- ;
			b1-- ;
			if ( a1 == orig_a1 ) break ;
		}
	}
	return prob ;
}


///////////////////////////////////////////////////////////////
// Chi-Square for a 2x2 table (approximation for fisher exact 
//                             2-tail if you have large values)
///////////////////////////////////////////////////////////////
//
// Calculates chi-square with an 2x2 table
// Arguments:  |   |   |
// 	       +--------
//             | a1| b1|
//             | a2| b2|
// See: http://faculty.vassar.edu/lowry/ch8pt2.html
///////////////////////////////////////////////////////////////


double chi_square_2x2( int a1, int a2, int b1, int b2 ) 
{
	int data[4] = { a1, a2, b1, b2 } ;
	int sums[4] ;
	double exp[4] ;
	int sum ;

	sum = data[0]+data[1]+data[2]+data[3] ;

	sums[0] = data[0]+data[1] ;
	sums[1] = data[2]+data[3] ;
	sums[2] = data[0]+data[2] ;
	sums[3] = data[1]+data[3] ;
	
	exp[0] = ((double)sums[0]*(double)sums[2])/(double)sum ;
	exp[1] = ((double)sums[0]*(double)sums[3])/(double)sum ;
	exp[2] = ((double)sums[1]*(double)sums[2])/(double)sum ;
	exp[3] = ((double)sums[1]*(double)sums[3])/(double)sum ;

	double result = 0. ;
	for ( int i=0 ; i < 4 ; ++i ) {
		double t = pow((double)(data[i])-exp[i], 2. ) / exp[i] ;
		result += t ;
	} 
	
	return  R::pchisq( result, 1., 0, 0 );
	return  0.03;
}	

// Calculates fisher exact if one value is below 10, otherwise chi_square
double fisher_chi2( int a1, int a2, int b1, int b2 ) {
	if ( a1 < 10 || a2 < 10 || b1 < 10 || b2 < 10 ) 
		return fishers_exact_2t( a1, a2, b1, b2 ) ;
	else 
		return chi_square_2x2( a1, a2, b1, b2 ) ;
}

#endif
