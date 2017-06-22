/*************************************************
*
* Copyright 2016 Boualem Boashash
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*
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
* Description:
*	
* Gateway Routine for zce.c
*************************************************/

#include <math.h>
#include "mex.h"
#include "matrix.h"

#include "arithm.h"		     /* routines for complex data types */
#include "tlocal.h"		     /* local function prototypes */

#include "zce.h"
#include "tfsa_c.h"

void mexFunction( int nlhs, MATRIX *plhs[], int nrhs, NEED_CONST MATRIX *prhs[])
{
   double *result;
   double *signal_r, *signal_i;
   int window_length;
   int signal_length;
   int M,N;

   /* check number of inputs */

   if( nrhs < 2 ) { /* zero or one input */

      tfsaErr( "zce", "Not enough input arguments" );
      winNTcheck( nlhs, plhs );
      return;
   }

   if( nrhs > 2 ) {	/* three or more inputs */

      tfsaErr( "zce", "Too many input arguments" );
      winNTcheck( nlhs, plhs );
      return;
   }

   /* Check first input -- can be complex */

   if( !mxIsNumeric( prhs[0] ) || !MXFULL( prhs[0] ) ) {
      tfsaErr( "zce", "Input must be a vector" );
      winNTcheck( nlhs, plhs );
      return;
   }

   /* Get input dimensions */

   M = ( int )mxGetM( prhs[0] );
   N = ( int )mxGetN( prhs[0] );

   /* Check input dimensions: A vector of unity length is considered
    * invalid; matrices are invalid. */

   if( (M==1 && N==1) || (M!=1 && N!=1) ) {
      tfsaErr( "zce", "Input must be a vector" );
      winNTcheck( nlhs, plhs );
      return;
   }

   signal_length = (int)(M>N?M:N);

   signal_r = mxGetPr( prhs[0] );
   if( mxIsComplex( prhs[0] ) )
     signal_i = mxGetPi( prhs[0]);
   else
     signal_i = NULL;


   /* Process second input */

   if( !GoodScalar( (MATRIX *)prhs[1] ) ) {
     tfsaErr( "zce", "Window length must be real scalar" );
     winNTcheck( nlhs, plhs );
     return;
   }

   window_length = (int)*(mxGetPr( prhs[1] )); /* get value */

   if( window_length < 1 ) {
     tfsaErr( "zce", "Window length must be greater than or equal to one" );
     winNTcheck( nlhs, plhs );
     return;
   }

   /* check window length */

   if( (window_length % 2) !=  0)
     window_length++;			    /* make it odd in length */

   if( signal_length < window_length ) {
      tfsaErr( "zce", "Window length must be less than or equal to signal length");
      winNTcheck( nlhs, plhs );
      return;
   }


   /* Check output */

   if( nlhs > 1 ) {
      tfsaErr( "zce", "Too many output arguments" );
      winNTcheck( nlhs, plhs );
      return;
   }

   plhs[0] = MXCREATEFULL( signal_length, 1, REAL );

   if( plhs[0] == NULL ) {
      tfsaErr( "zce", "Memory allocation failed" );
      winNTcheck( nlhs, plhs );
      return;
   }
   result = mxGetPr( plhs[0] );

   /* Zero-Crossing Estimator */

   zce( signal_r, signal_i, signal_length, result, window_length );

}
