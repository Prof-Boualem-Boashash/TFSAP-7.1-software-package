/*************************************************
* Author:                 Boualem Boashash         (boualem.boashash@gmail.com)
* Maintainer since 2015:  Samir Ouelha  			(samir_ouelha@hotmail.fr)
*
* The following 2 references should be cited whenever this script is used:
* [1] B. Boashash, Samir Ouelha, Designing time-frequency  and time-scale
* features for efficient classification of non stationary signals: a tutorial
* review with performance comparison, Digital Signal Processing, 2017.
* [2] B. Boashash, Samir Ouelha, Efficient Software Matlab Package for the calculation
* of Time-Frequency Distributions related Time-Scale methods and the extraction of
* signal characteristics, SoftwareX, 2017.
* In addition, the following 3rd reference is relevant for use of the executable
* [3] B. Boashash (ed.), Time-Frequency Signal Analysis and Processing, 2nd Edition
* (London: Elsevier / Academic Press, December 2015); ISBN 978-0-12-398499-9. Book website:
* http://booksite.elsevier.com/9780123984999/
*
*
* arithm support function
* 
*************************************************/

#include   <math.h>

#include   "arithm.h"

complex	cmplx (real_part, imaginary_part)
 double	real_part, imaginary_part;
/*--------------------------------------
 * Assign value	to a complex number
 * Input: two real numbers;
 * Output: Complex number
 *--------------------------------------*/
{
  complex x;

    x.re = real_part;
    x.im = imaginary_part;
    return (x);
}


complex	add (x,	y)
 complex x, y;
/*---------------------------------------
 * Sum of two complex numbers
 *_______________________________________*/
{
    complex z;

    z.re = x.re	+ y.re;
    z.im = x.im	+ y.im;
    return (z);
}


complex	sub (x,	y)
 complex x, y;
/*---------------------------------------
 * Substuct two	complex	numbers
 *--------------------------------------*/
{
    complex z;

    z.re = x.re	- y.re;
    z.im = x.im	- y.im;
    return (z);
}


complex	multpl (x, y)
 complex x,y;
/*---------------------------------------
 * Multiply two	complex	numbers
 *---------------------------------------*/
{
   complex z;
   z.re	= x.re * y.re -	x.im * y.im;
   z.im	= x.re * y.im +	x.im * y.re;
   return (z);
}


complex divide( complex x, complex y )
{
/*---------------------------------------
 * Divide two complex numbers x/y
 *---------------------------------------*/
  double denom = 0;
  complex nom;
  double E = 1e-14;

  /* If denominator or nominator is nearly zero then let it equal zero so as */
  /* to avoid inf. and divide by zero errors... */
  if( ( fabs(y.re) < E ) && ( fabs(y.re) > 0 ) )
    y.re = 0.0;
  if( ( fabs(y.im) < E ) && ( fabs(y.im) > 0 ) )
    y.im = 0.0;
  if( ( fabs(x.re) < E ) && ( fabs(x.re) > 0 ) )
    x.re = 0;
  if( ( fabs(x.im) < E ) && ( fabs(x.im) > 0 ) )
    x.im = 0;

  
  denom = (( y.re * y.re ) + ( y.im * y.im ));
  /* Should throw error here.. */
  if( denom == 0 )
    return( cmplx( (double)0, (double)0 ) );

  nom = multpl( x, conj( y ) );

  return( cmplx( nom.re / denom , nom.im / denom ) );
}

  



complex	conj (x)
  complex x;
/*---------------------------------------------
 * Complex-conjugate of	x
 *--------------------------------------------*/
{
 complex z;

 z.re =	x.re;
 z.im =	-x.im;
 return	(z);
}



double module(x)
 complex x;
/*-------------------------------------
 * Absolute Value of complex number
 *-------------------------------------*/
{
 double	z;

  z = x.re * x.re + x.im * x.im	;
  return (sqrt(z));

}


double arg (x)
  complex x;
{
  double theta;

  theta	= atan (x.im/x.re);
  return (theta);

}


complex	cmult (x,c)
  complex x;
  double c;
/*------------------------------------------
 * multiply complex number by real constant
 *------------------------------------------*/
{

  complex z;
  z.re = c * x.re;
  z.im = c * x.im;
  return (z);
}




