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
* tlocal.c: TFSA local	functions
* 
*************************************************/
#include <math.h>
#include "mex.h"
#include <stdio.h>

#include "arithm.h"		     /* routines for complex data */
#include "tlocal.h"		     /* local function prototypes */
#include "tfsa_c.h"

void tfsaErr( char *name, char *string )
{
  mexPrintf( "??? Error using ==> %s%c\n", name, 7 );
  mexPrintf( "%s.\n\n", string );

}

void tfsaWarning( char *name, char *string )
{
  mexPrintf( "??? Warning ==> %s\n", name );
  mexPrintf( "%s.\n\n", string );

}

/* void winNTcheck( int nlhs, Matrix *plhs[] )
 *
 * winNTcheck should be called every time the gateway returns with an
 * error condition.  An error exists with Mathwork's mex interface
 * under Windows NT, causing a general protection fault (GPF) error
 * whenever a function is called with no input arguments but with an
 * output argument(s). As a solution to this problem, all output
 * matrices are allocated a null matrix.
 *
 * NB: Since this problem does not exist under UNIX or Windows 3.1,
 * the constant WINDOWS_NT is defined to unnecessarily producing
 * output.  */

void winNTcheck( int nlhs, MATRIX *plhs[] ) {

#ifdef WINDOWS_NT

  int i;

  for( i=0; i < nlhs; i++ )
    plhs[i] = MXCREATEFULL( 0, 0, REAL );   /* NULL array */

#endif

}


/* Type complex conversion functions */

void ToComplex( complex *vector, double *re, double *im, int N )
{

  int i;

  for( i=0; i < N; i++ )  {
    vector[i].re = re[i];
    vector[i].im = im[i];
  }

 return;

}

void RealToComplex( complex *vector, double *re, int N )
{

  int i;

  for( i=0; i < N; i++ )
    vector[i].re = re[i];


 return;

}


void FromComplex( complex *vector, double *re, double *im, int N )
{

  int i;

  for( i=0; i < N; i++ )  {
    re[i] = vector[i].re;
    im[i] = vector[i].im;
  }

 return;

}


/* GoodScalar
 *
 * Check to see if given pointer is a non-negative scalar. Returns:
 *
 * 0 on fail (a scalar, complex scalar, matrix, or non-numeric input)
 *
 * 1 on pass (otherwise ok)
 */

int GoodScalar( MATRIX *ptr )
{
  int M, N;

  if( !mxIsNumeric( ptr ) || !MXFULL( ptr ) )
    return( 0 );

  /* Check dimensions*/

  M = ( int )mxGetM( ptr );
  N = ( int )mxGetN( ptr );

  /*  printf( "M=%d, N=%d\n", M, N ); */

  if( (M==1) && (N==1) )
    return 1;

  return 0;
}

