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
*  Gateway Routine for rihaczek.c
*************************************************/

#include <math.h>
#include <string.h>

#include "mex.h"
#include "matrix.h"

#include "arithm.h"		     /* routines for complex data */
#include "tlocal.h"		     /* local function prototypes */

#include "rihaczek.h"
#include "tfsa_c.h"
#include "window.h"

void mexFunction( int nlhs, MATRIX *plhs[], int nrhs, NEED_CONST MATRIX *prhs[])
{
  double *result_re, *result_im;
  double *signal_real, *signal_imag;
  int signal_length, signal_r2, signal_order;
  int M,N,N_min;
  int fft_length, numslices, time_res;
  int type_of_dist, n;
  int window_length, window_type;
  char *p;


  /* basic input--output number check */
  
  if( nrhs < 1 ) {
    tfsaErr( "rihaczek", "Not enough input arguments" );
    winNTcheck( nlhs, plhs );
    return;
  }
  
  if( nrhs > 7 ) {
    tfsaErr( "spramn", "Too many input arguments" );
    winNTcheck( nlhs, plhs );
    return;
  }
  
  if( nlhs > 1 ) {
    tfsaErr( "rihaczek", "Too many output arguments" );
    winNTcheck( nlhs, plhs );
    return;
  }

  
  
  /* Check first input */
  
  if( !mxIsNumeric( prhs[0] ) || !MXFULL( prhs[0] ) ) {
    tfsaErr( "rihaczek", "Input must be a vector" );
    winNTcheck( nlhs, plhs );
    return;
  }
  
  
  /* Get input matrix dimensions */
  
  M = ( int )mxGetM( prhs[0] );
  N = ( int )mxGetN( prhs[0] );
  
  /* Check input dimensions: A vector of unity length is considered
   * invalid; matrices are invalid. */
  
  if( (M==1 && N==1) || (M!=1 && N!=1) ) {
    tfsaErr( "rihaczek", "Input must be a vector" );
    winNTcheck( nlhs, plhs );
    return;
  }
  
  signal_length = (int)(M>N?M:N);

  /* Check second input */
  if( nrhs > 1 ) {
    if( !GoodScalar( (MATRIX *)prhs[1] ) ) {
      tfsaErr( "rihaczek", "Time resolution must be a scalar" );
      winNTcheck( nlhs, plhs );
      return;
    }
    else
      time_res = (int)*(mxGetPr( prhs[1] ));
  } else {
    time_res = 1;
  }

  /* Parameter value check */
  if( time_res < 1 ) {
    tfsaErr( "rihaczek", "Time resolution length must be greater than zero" );
    winNTcheck( nlhs, plhs );
    return;
  }

  
  /* Check third input */
  if( nrhs > 2 ){
    if( !GoodScalar( (MATRIX *)prhs[2] ) ) {
      tfsaErr( "rihaczek", "FFT length must be a scalar" );
      winNTcheck( nlhs, plhs );
      return;
    }
    else
      fft_length = (int)*(mxGetPr( prhs[2] ) );  /* Get value */
  } else {
    fft_length = signal_length;
  }
    
  /* Parameter value check */
  if( fft_length < 1 ) {
    tfsaErr( "rihaczek", "FFT length must be greater than zero" );
    winNTcheck( nlhs, plhs );
    return;
  }
  
  
  /* calculate radix-2 value equal or above fft_length */
  signal_order = 0;
  signal_r2 = 1;
  while( signal_r2 < fft_length){
      signal_order++;
      signal_r2 <<= 1;
  }

  /* Check fourth input */
  if( nrhs > 3 ){
    if( !GoodScalar( (MATRIX *)prhs[3] ) ) {
      tfsaErr( "rihaczek", "Rihaczek or Levin option must be a scalar" );
      winNTcheck( nlhs, plhs );
      return;
    }
    else
      type_of_dist = (int)*(mxGetPr( prhs[3] ) );  /* Get value */
  } else {
    type_of_dist = 0;
  }

  
  /* Parameter value check */
  if( (type_of_dist != 1) && (type_of_dist != 0) ) {
    tfsaErr( "rihaczek", "1 or 0 required for Rihaczek = 0 and Levin = 1" );
    winNTcheck( nlhs, plhs );
    return;
  }


  signal_real = mxGetPr( prhs[0]);
  if (mxIsComplex(prhs[0]))
    signal_imag = mxGetPi( prhs[0]);
  else
    signal_imag = NULL;
  

  /* Find out which is smaller, FFT length or signal length */
  N_min = (int) (signal_r2 < signal_length ? signal_r2 : signal_length );


  /* Calculate number of times for computing FFT. This
   * is equal to the number of estimated frequency points
   */
  numslices = (int) ceil((double)N_min/time_res);


  /* If Levin distribution then just take the real part */  
  if( type_of_dist == 0 )
    plhs[0] = MXCREATEFULL ((signal_r2/2)+1, numslices, COMPLEX);
  else
    plhs[0] = MXCREATEFULL ((signal_r2/2)+1, numslices, REAL);

  if (plhs[0] == NULL)
    {
      tfsaErr("rihaczek", "Memory allocation failed" );
      winNTcheck( nlhs, plhs );
      return;
    }
  result_re = mxGetPr( plhs[0]);

  if( type_of_dist == 0 )
    result_im = mxGetPi( plhs[0]);
  else
    result_im = NULL;



  /* Check fifth input */
  if( nrhs > 4 ){
    if( !GoodScalar( (MATRIX *)prhs[4] ) ) {
      tfsaErr( "rihaczek", "Smoothing window length must be a scalar" );
      winNTcheck( nlhs, plhs );
      return;
    }
    else
      window_length = (int)*(mxGetPr( prhs[4] ) );  /* Get value */
  
    /* Parameter value check */
    if( window_length < 1 ) {
      tfsaErr( "rihaczek", "Window length must be greater than zero" );
      winNTcheck( nlhs, plhs );
      return;
    }
  
    /* Check window length against signal length */
    if( window_length > signal_length ) {
      tfsaWarning( "rihaczek", "Window length has been truncated to signal length" );
      window_length = (int)signal_length;
    }
  
    /* Check for odd window length.. */
    if( fmod( window_length, 2) == 0 ){
      window_length = (int)(window_length - 1);
      tfsaWarning( "rihaczek", "Window length is not odd, truncating to odd value." );
    } 

  } else {
    window_length = 0;
  }

  
  /* Check sixth input*/
  if( nrhs > 5 ){
    if( MXSTRING( prhs[5] ) ) {
      
      
      n = (int)mxGetN( prhs[5] ) + 1;  /* string length + 1 for NULL */
      p = mxCalloc( n, sizeof( char ) );
      if( p == NULL ) {
	tfsaErr( "rihaczek", "Internal memory allocation failure" );
	winNTcheck( nlhs, plhs );
	return;
      }
      
      if( mxGetString( prhs[5], p, n ) ) {
	tfsaErr( "rihaczek", "Could not get smoothing window type" );
	winNTcheck( nlhs, plhs );
	return;
      }
      
      if( !strcmp( p, "rect" ) || !strcmp( p, "RECT" ) )
	window_type = RECT;
      else if( !strcmp( p, "hann" ) || !strcmp( p, "HANN" ) )
	window_type = HANN;
      else if( !strcmp( p, "hamm" ) || !strcmp( p, "HAMM" ) )
	window_type = HAMM;
      else if( !strcmp( p, "bart" ) || !strcmp( p, "BART" ) )
	window_type = BART;
      else {
	tfsaErr( "rihaczek", "Unknown smoothing window type" );
	winNTcheck( nlhs, plhs );
	return;
      }
    } else {
      tfsaErr( "rihaczek",  "Smoothing window type must be a string" );
      winNTcheck( nlhs, plhs );
      return;
    }
  } else {
    window_type = 0;
  }
  if( nrhs == 5 ){
    tfsaErr( "rihaczek",  "Need to specify window type" );
    winNTcheck( nlhs, plhs );
    return;
  }


  rihaczek( signal_real, signal_imag, signal_length, result_re, result_im,  
  	    numslices, time_res, signal_r2, signal_order, type_of_dist,  
  	    window_length, window_type );   

}


