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
* Description:
*	
* gateway routine for pdeIF.c.
*************************************************/
#include <math.h>
#include "mex.h"
#include <matrix.h>

#include "arithm.h"		     /* routines for complex data */
#include "tlocal.h"		     /* local function prototypes */

#include "pdeIF.h"
#include "tfsa_c.h"

void mexFunction( int nlhs, MATRIX *plhs[], int nrhs, NEED_CONST MATRIX *prhs[])
{
   double *result;
   double *signal_r, *signal_i;
   int window_length;
   int signal_length;
   int M,N;
   int estimator_order, Kay_smoothing;


   /* basic input--output number check */

   if( nrhs < 2 ) {
      tfsaErr( "pde", "Not enough input arguments" );
      winNTcheck( nlhs, plhs );
      return;
   }

   if( nrhs > 3 ) {
      tfsaErr( "pde", "Too many input arguments" );
      winNTcheck( nlhs, plhs );
      return;
   }

   if( nlhs > 1 ) {
      tfsaErr( "pde", "Too many output arguments" );
      winNTcheck( nlhs, plhs );
      return;
   }

   /* Check first input */


   if( !mxIsNumeric( prhs[0] ) || !MXFULL( prhs[0] ) ) {
      tfsaErr( "pde", "Input must be a vector" );
      winNTcheck( nlhs, plhs );
      return;
   }

   /* Get input matrix dimensions */

   M = ( int )mxGetM( prhs[0] );
   N = ( int )mxGetN( prhs[0] );

   /* Check input dimensions: A vector of unity length is considered
    * invalid; matrices are invalid. */

   if( (M==1 && N==1) || (M!=1 && N!=1) ) {
      tfsaErr( "pde", "Input must be a vector" );
      winNTcheck( nlhs, plhs );
      return;
   }

   /* Choose the larger row or column dimension */

   signal_length = (int)( M>N ? M:N );


   /* Check second input */
   if( !GoodScalar( (MATRIX *)prhs[1] ) ) {
     tfsaErr( "pde", "Estimator order must be a scalar" );
     winNTcheck( nlhs, plhs );
     return;
   }
   else
     estimator_order = (int)*(mxGetPr( prhs[1]));


   /* Parameter value check */
   if( estimator_order < 1 ) {
     tfsaErr( "pde", "Estimator order must be greater than zero" );
     winNTcheck( nlhs, plhs );
     return;
   }

  if( (estimator_order != 1) && (estimator_order != 2) &&
      (estimator_order != 4) && (estimator_order != 6) ) {
     tfsaErr( "pde", "Estimator order must 1, 2, 4, or 6" );
     winNTcheck( nlhs, plhs );
     return;
   }
 
   if (nrhs == 3) { 
     /* Check second input */

     if( !GoodScalar( (MATRIX *)prhs[2] ) ) {
       tfsaErr( "pde", "Smoothing window length must be a scalar" );
       winNTcheck( nlhs, plhs );
       return;
     }
     else
       window_length = (int)*(mxGetPr( prhs[2] ) );  /* Get value */

     /* Parameter value check */
     if( window_length < 1 ) {
        tfsaErr( "pde", "Window length must be greater than zero" );
        winNTcheck( nlhs, plhs );
        return;
     }

     /* Check window length against signal length */

     if( signal_length < (2*window_length) ) {
       tfsaWarning( "pde", "Window length has been truncated to half signal length" );
       window_length = (int)(signal_length/2);
     }
     Kay_smoothing = 1;
   } else {
     window_length = 0;
     Kay_smoothing = 0; 
   }

  
   /* Dereference input */

   signal_r = mxGetPr( prhs[0] );
   if( mxIsComplex( prhs[0] ) )
     signal_i = mxGetPi( prhs[0] );
   else
     signal_i = NULL;


   /* Create output vector */

   plhs[0] = MXCREATEFULL(signal_length, 1, REAL);

   if( plhs[0] == NULL ) {
      tfsaErr( "pde", "Internal memory allocation failure" );
      winNTcheck( nlhs, plhs );
      return;
   }

   result = mxGetPr( plhs[0] );

   /* Phase Difference Estimator and smoothed Kay estimator*/

   if( pdeIF(signal_r, signal_i, signal_length, result, window_length,
	   estimator_order, Kay_smoothing) )
     winNTcheck( nlhs, plhs );  /* An error occured in pdeIF.c */

}
