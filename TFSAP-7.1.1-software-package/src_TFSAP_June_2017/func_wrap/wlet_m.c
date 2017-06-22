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
* gateway routine for wlet.c.
*************************************************/
#include <math.h>
#include "mex.h"
#include "matrix.h"

#include "arithm.h"		     /*	routines for complex data types */
#include "tlocal.h"		     /* local function prototypes */

#include "wlet.h"
#include "window.h"
#include "tfsa_c.h"

void mexFunction( int nlhs, MATRIX *plhs[], int nrhs, NEED_CONST MATRIX *prhs[])
{
   double *result, *result2;
   double *signal_r;
   int signal_length;
   int M,N;
   int num_coeff;
   int stages;
   int output_type;
   int signal_order, signal_r2;
   int direction;


   /* check input arguments */

   if( nrhs == 0 ) {  /* No input */

      tfsaErr( "wlet", "Not enough input arguments" );
      winNTcheck( nlhs, plhs );
      return;
   }

   if( nrhs > 4 ) {	/* Five or more inputs */

      tfsaErr( "wlet", "Too many input arguments" );
      winNTcheck( nlhs, plhs );
      return;
   }

   /* 1 to 4 inputs given */

   /* Check first input */

   if( !mxIsNumeric( prhs[0] ) || !MXFULL( prhs[0] ) ||
      mxIsComplex( prhs[0] ) ) {
     tfsaErr( "wlet", "Input must be a real vector" );
     winNTcheck( nlhs, plhs );
     return;
   }

   /* Get input dimensions */

   M = ( int )mxGetM( prhs[0] );
   N = ( int )mxGetN( prhs[0] );

   /* Check input dimensions: A vector of unity length is considered
    * invalid; matrices are invalid. */

   if( (M==1 && N==1) || (M!=1 && N!=1) ) {
      tfsaErr( "wlet", "Input must be a vector" );
      winNTcheck( nlhs, plhs );
      return;
   }

   signal_length = (int)(M>N?M:N);


   /* calculate radix-2 value equal or above window_length */
   signal_order = 0;
   signal_r2 = 1;

   while( signal_r2 < signal_length) {
      signal_order++;
      signal_r2 <<= 1;
   }

   if( signal_length != signal_r2 ) {
      tfsaErr( "wlet", "Signal length must be a radix-2 number" );
      winNTcheck( nlhs, plhs );
      return;
   }

   /* Dereference first input */

   signal_r = mxGetPr( prhs[0] );



   /* Process the second input if it exists, otherwise assume default */

   if( nrhs > 1 ) {
      if( !GoodScalar( (MATRIX *)prhs[1] ) ) {
         tfsaErr( "wlet", "Output type must be 1 or 2" );
	 winNTcheck( nlhs, plhs );
	 return;
      }

      output_type = (int)*(mxGetPr( prhs[1] ) );  /* Get value */

      if( (output_type != ONE_D) && (output_type != TWO_D) ) {
	 tfsaErr( "wlet", "Output type must be 1 or 2" );
	 winNTcheck( nlhs, plhs );
	 return;
      }
   }
   else
     output_type = ONE_D;		   /* Default */

   /* Process the third input if it exists, otherwise assume default */

   if( nrhs > 2 ) {
      if( !GoodScalar( (MATRIX *)prhs[2] ) ) {
	 tfsaErr( "wlet", "The number of coefficients must be 4, 12, or 20" );
	 winNTcheck( nlhs, plhs );
	 return;
      }

      num_coeff = (int) *(mxGetPr( prhs[2] )); /* Get value */

      if( (num_coeff != 4) && (num_coeff != 12) && (num_coeff != 20)) {
	 tfsaErr( "wlet", "The number of coefficients must be 4, 12, or 20" );
	 winNTcheck( nlhs, plhs );
	 return;
      }
   }
   else
     num_coeff = 20;			/* Default */

   /* Process the fourth input if it exists, otherwise assume default */

   if( nrhs > 3 ) {
      if( !GoodScalar( (MATRIX *)prhs[3] ) ) {
	 tfsaErr( "wlet", "The direction flag must be -1 or 1" );
	 winNTcheck( nlhs, plhs );
	 return;
      }

      direction = (int) *(mxGetPr( prhs[3] ));

      if( (direction != -1 ) && (direction != 1 ) ) {
	 tfsaErr( "wlet", "The direction flag must be -1 or 1" );
	 winNTcheck( nlhs, plhs );
	 return;
      }
   }
   else
     direction = 1; 	      /* Default = forward */


   if( (direction == -1) && (output_type == TWO_D) ) {
      tfsaErr( "wlet", "Two Dimensional output is incompatible with reverse transform");
      winNTcheck( nlhs, plhs );
      return;
   }

   /* Check output arguments.  All input parameters are ok at this stage. */


   if( nlhs > 1 ) {
      tfsaErr( "wvd", "Too many output arguments" );
      winNTcheck( nlhs, plhs );
      return;
   }


   if( output_type == ONE_D ) {	   /* 1-D transform */

      plhs[0] = MXCREATEFULL( signal_length, 1, REAL );

      if( plhs[0] == NULL ) {
	 tfsaErr( "wvd", "Memory allocation failed" );
	 winNTcheck( nlhs, plhs );
	 return;
      }

      result = mxGetPr( plhs[0] );

      wave( signal_r, result, signal_length, direction, num_coeff );
   }
   else {	    /* 2-D transform */

      result = (double *)mxCalloc( signal_length, sizeof(double) );  /* temp */
      if( result == NULL ) {
	 tfsaErr( "wvd", "Memory allocation failed" );
	 winNTcheck( nlhs, plhs );
	 return;
      }

      /* Perform the transform */
      /* Note that only forward transform makes sense here */

      stages = wave( signal_r, result, signal_length, 1, num_coeff);
      if( stages < 0 ) {
	 tfsaErr( "wvd", "Failed due to stages < 0" ); /* shouldn't happen */
         winNTcheck( nlhs, plhs );
	 return;
      }

      plhs[0] = MXCREATEFULL( signal_length/2, signal_length, REAL);

      if( plhs[0] == NULL ) {
         tfsaErr( "wvd", "Memory allocation failure" );
         winNTcheck( nlhs, plhs );
	 return;
      }

      result2 = mxGetPr( plhs[0] );   /* actual result */

      /* Convert result into a 2D matrix */

      form_ts( result, result2, signal_length );
   }

}

