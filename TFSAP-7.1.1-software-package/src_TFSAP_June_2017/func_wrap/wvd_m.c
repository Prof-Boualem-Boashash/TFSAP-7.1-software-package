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
* Gateway function for WVD
*************************************************/

#include <math.h>
#include "mex.h"
#include "matrix.h"

#include "arithm.h"		     /* routines for complex data */
#include "tlocal.h"		     /* local function prototypes */

#include "wvd.h"
#include "tfsa_c.h"

void mexFunction( int nlhs, MATRIX *plhs[], int nrhs, NEED_CONST MATRIX *prhs[])
{
  double *result;		       /* pointer to results matrix data */
  double *signal_r, *signal_i;
  int window_length, time_res;
  int signal_length;
  int M,N;
  int fft_length;
  int num_slices, window_order, window_r2;


  /* basic input--output number check */

   if( nrhs < 1 ) {
     tfsaErr( "wvd", "Not enough input arguments" );
     winNTcheck( nlhs, plhs );
     return;
   }

  if( nrhs > 4 ) {
    tfsaErr( "wvd", "Too many input arguments" );
    winNTcheck( nlhs, plhs );
      return;
  }

  if( nlhs > 1 ) {
    tfsaErr( "wvd", "Too many output arguments" );
    winNTcheck( nlhs, plhs );
    return;
  }

  /* Check first input */

  if( !mxIsNumeric( prhs[0] ) || !MXFULL( prhs[0] ) ) {
    tfsaErr( "wvd", "Input must be a vector" );
    winNTcheck( nlhs, plhs );
    return;
   }


  /* Get input matrix dimensions */

  M = ( int )mxGetM( prhs[0] );
  N = ( int )mxGetN( prhs[0] );

  /* Check input dimensions: A vector of unity length is considered
   * invalid; matrices are invalid. */

  if( (M==1 && N==1) || (M!=1 && N!=1) ) {
    tfsaErr( "wvd", "Input must be a vector" );
    winNTcheck( nlhs, plhs );
    return;
  }

  signal_length = (unsigned)(M>N?M:N);

  /* Do defaults */

  switch( nrhs ) {

  case( 1 ):  /* Install all defaults */

      window_length = DEF_WVD_WIN_LEN;
      fft_length = window_length;
      time_res = 1;
      break;

    case( 2 ):  /* Get window length */

      if( !GoodScalar( (MATRIX *)prhs[1] ) ) {
	tfsaErr( "wvd", "Window length must be a scalar" );
	winNTcheck( nlhs, plhs );
	return;
      }
      else
	window_length = (int)*(mxGetPr( prhs[1] ) );  /* Get value */

      /* Parameter value check */
      if(  window_length < 1 ) {
	tfsaErr( "wvd", "Window length must be greater than zero" );
	winNTcheck( nlhs, plhs );
	return;
      }

      /* Check window length against signal length */

      if( window_length > signal_length )  {
	 tfsaWarning( "wvd", "Window length has been truncated to signal length" );
         if ((signal_length%2) == 1) {
	    window_length = (int) signal_length;
         } else {
            window_length = (int) signal_length-1;
         }
      }

      if( (window_length%2) == 0) {
	tfsaErr("wvd", "Window length must be odd" );
	winNTcheck( nlhs, plhs );
        return;
      }

      fft_length = window_length+1;
      time_res = 1;

      break;

    case( 3 ):

      if( !GoodScalar( (MATRIX *)prhs[1] ) ) {
	tfsaErr( "wvd", "Window length must be a scalar" );
	winNTcheck( nlhs, plhs );
	return;
      }
      else
	window_length = (int)*(mxGetPr( prhs[1] ) );  /* Get value */

      /* Parameter value check */
      if(  window_length < 1 ) {
	tfsaErr( "wvd", "Window length must be greater than zero" );
	winNTcheck( nlhs, plhs );
	 return;
      }

      /* Check window length against signal length */

      if( window_length > signal_length )  {
 	 tfsaWarning( "wvd", "Window length has been truncated to signal length" );
         if ((signal_length%2) == 1) {
	    window_length = (int) signal_length;
         } else {
            window_length = (int) signal_length-1;
         }
      }

      if( (window_length%2) == 0) {
	tfsaErr("wvd", "Window length must be odd" );
	winNTcheck( nlhs, plhs );
	return;
      }


      /* Check third input */
      if( !GoodScalar( (MATRIX *)prhs[2] ) ) {
	tfsaErr( "wvd", "Time resolution must be a scalar" );
	winNTcheck( nlhs, plhs );
	return;
      }
      else
	time_res = (int)*(mxGetPr( prhs[2]));

      /* Parameter value check */
      if( time_res < 1 ) {
	tfsaErr( "wvd", "Time resolution must be greater than zero" );
         winNTcheck( nlhs, plhs );
	return;
      }

      if (time_res > signal_length) {
         tfsaErr( "wvd", "Time resolution must be no greater than signal length" );
         winNTcheck( nlhs, plhs );
         return;
      }

      fft_length = window_length;

      break;

    case( 4 ):


      if( !GoodScalar( (MATRIX *)prhs[1] ) ) {
	tfsaErr( "wvd", "Window length must be a scalar" );
	winNTcheck( nlhs, plhs );
	return;
      }
      else
	window_length = (int)*(mxGetPr( prhs[1] ) );  /* Get value */

      /* Parameter value check */
      if(  window_length < 1 ) {
	tfsaErr( "wvd", "Window length must be greater than zero" );
	winNTcheck( nlhs, plhs );
	return;
      }

      /* Check window length against signal length */

      if( window_length > signal_length )  {
	 tfsaWarning( "wvd", "Window length has been truncated to signal length" );
         if ((signal_length%2) == 1) {
	    window_length = (int) signal_length;
         } else {
            window_length = (int) signal_length-1;
         }
      }

      if( (window_length%2) == 0) {
	tfsaErr("wvd", "Window length must be odd" );
	winNTcheck( nlhs, plhs );
	return;
      }


      /* Check third input */
      if( !GoodScalar( (MATRIX *)prhs[2] ) ) {
	tfsaErr( "wvd", "Time resolution must be a scalar" );
	winNTcheck( nlhs, plhs );
	return;
      }
      else
	time_res = (int)*(mxGetPr( prhs[2]));

      /* Parameter value check */
      if( time_res < 1 ) {
	tfsaErr( "wvd", "Time resolution must be greater than zero" );
	winNTcheck( nlhs, plhs );
	return;
      }

      if (time_res > signal_length) {
        tfsaErr( "wvd", "Time resolution must be no greater than signal length" );
        winNTcheck( nlhs, plhs );
        return;
      }

      /* Check fourth input */

      if( !GoodScalar( (MATRIX *)prhs[3] ) ) {
	tfsaErr( "wvd", "FFT length must be a scalar" );
	winNTcheck( nlhs, plhs );
	return;
      }
      else
	fft_length = (int)*(mxGetPr( prhs[3]));


      if( fft_length < 1 ) {
	tfsaErr( "wvd", "FFT length must be greater than zero" );
	winNTcheck( nlhs, plhs );
	return;
      }
      if (fft_length < window_length)
	fft_length = window_length;

      break;
    }


   num_slices = (int) ceil((double)signal_length/time_res);

   /*  mexPrintf("num_slices is:  %d",num_slices); */

   /* calculate radix-2 value equal or above window_length */
   window_order = 0;
   window_r2 = 1;
   while( window_r2 < fft_length) {
      window_order++;
      window_r2 <<= 1;
   }

   /* Allocate memory for the output variable result
    *
    * NB: The WVD is always real, whether the input is real or complex.
    *
    */

   plhs[0] = MXCREATEFULL( window_r2, num_slices, REAL );
   if (plhs[0] == NULL) {
      tfsaErr("wvd", "Memory allocation failed" );
      winNTcheck( nlhs, plhs );
      return;
   }
   result = mxGetPr( plhs[0] );


   /* Dereference input */

   signal_r = mxGetPr( prhs[0] );
   if( mxIsComplex(prhs[0]) )
     signal_i = mxGetPi( prhs[0] );
   else
     signal_i = NULL;


  wvd( signal_r, signal_i, signal_length, result, num_slices,
      time_res, window_length, window_r2, window_order);
}

