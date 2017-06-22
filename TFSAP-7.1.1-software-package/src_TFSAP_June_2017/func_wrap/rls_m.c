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
* Gateway Routine for rls.c
*/

#include <math.h>

#include "mex.h"
#include "matrix.h"

#include "arithm.h"		     /* routines for complex data types */
#include "tlocal.h"		     /* local function prototypes */

#include "rls.h"
#include "tfsa_c.h"

void mexFunction( int nlhs, MATRIX *plhs[], int nrhs, NEED_CONST MATRIX *prhs[])
{
   double *result;
   double *signal_r, *signal_i;
   int signal_length;
   int M,N;
   float alpha;

   /* check input arguments */

   switch( nrhs ) {

    case( 0 ):  /* No input */
    case( 1 ):  /* One input */

      tfsaErr("rls", "Not enough input arguments" );
      winNTcheck( nlhs, plhs );
      return;
      break;

    case( 2 ):  /* Two inputs */

      /* Check first input  -- can be complex */

      if( !mxIsNumeric( prhs[0] ) || !MXFULL( prhs[0] ) ) {
	 tfsaErr( "rls", "Input must be a vector" );
	 winNTcheck( nlhs, plhs );
	 return;
      }

      /* Get input dimensions */

      M = ( int )mxGetM( prhs[0] );
      N = ( int )mxGetN( prhs[0] );

      /* Check input dimensions: A vector of unity length is considered
       * invalid; matrices are invalid. */

      if( (M==1 && N==1) || (M!=1 && N!=1) ) {
	 tfsaErr( "rls", "Input must be a vector" );
	 winNTcheck( nlhs, plhs );
	 return;
      }

      /* Choose the larger row or column dimension */

      signal_length = (int)(M>N?M:N);

      /* Dereference first input */

      signal_r = mxGetPr( prhs[0] );
      if( mxIsComplex( prhs[0] ) )
	signal_i = mxGetPi( prhs[0] );
      else
	signal_i = NULL;

      /* check second input */

      if( !GoodScalar( (MATRIX *)prhs[1] ) ) {
	 tfsaErr( "rls", "Forgetting factor 'alpha' must be a real scalar" );
	 winNTcheck( nlhs, plhs );
	 return;
      }

      alpha = (float)*(mxGetPr( prhs[1] ) );	/* Get value */

      if( alpha < 0.0 ) {
	 tfsaErr( "rls", "Forgetting factor 'alpha' must be a non-negative scalar" );
	 winNTcheck( nlhs, plhs );
	 return;
      }
      break;

  default:   /* More than two inputs */

      tfsaErr( "rls", "Too many input arguments" );
      winNTcheck( nlhs, plhs );
      return;
      break;

   }  /* end of switch( nrhs ) */



   /* Check output arguments.  All input parameters are ok at this stage. */

   switch( nlhs ) {

    case( 0 ):
    case( 1 ):

      /* allocate output row vector, and check memory allocation */

      plhs[0] = MXCREATEFULL( signal_length, 1, REAL );

      if( plhs[0] == NULL ) {
	 tfsaErr( "rls", "Internal memory allocation failure" );
	 winNTcheck( nlhs, plhs );
	 return;
      }

      result = mxGetPr( plhs[0] );

      break;

    default:  /* more than one output */

      tfsaErr( "rls", "Too many output arguments" );
      winNTcheck( nlhs, plhs );
      return;
      break;
   }


   /* Recursive Least Square Adaptive Estimator */

   rls( signal_r, signal_i, signal_length, result, alpha );

}
