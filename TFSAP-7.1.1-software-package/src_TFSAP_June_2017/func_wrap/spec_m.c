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
* gateway routine for spec.c.
*************************************************/

#include <math.h>
#include <string.h>

#include "mex.h"
#include "matrix.h"

#include "arithm.h"		     /* routines for complex data */
#include "tlocal.h"		     /* local function prototypes */

#include "spec.h"
#include "tfsa_c.h"
#include "window.h"

void mexFunction( int nlhs, MATRIX *plhs[], int nrhs, NEED_CONST MATRIX *prhs[])
{
  double *tfd_re, *tfd_im;
  double *signal_real, *signal_imag;
  int window_length, time_res;
  int signal_length;
  int M,N,N_min;
  int fft_length;
  int numslices, window_order, window_r2;
  int window_type;
  unsigned stft_or_spec;

  int n;
  char *p;

  /* basic input--output number check */
  
  if( nrhs < 4 ) {
    tfsaErr( "spec", "Not enough input arguments" );
    winNTcheck( nlhs, plhs );
    return;
  }
  
  if( nrhs > 6 ) {
    tfsaErr( "spec", "Too many input arguments" );
    winNTcheck( nlhs, plhs );
    return;
  }
  
  if( nlhs > 1 ) {
    tfsaErr( "spec", "Too many output arguments" );
    winNTcheck( nlhs, plhs );
    return;
  }

  
  
  /* Check first input */
  
  if( !mxIsNumeric( prhs[0] ) || !MXFULL( prhs[0] ) ) {
    tfsaErr( "spec", "Input must be a vector" );
    winNTcheck( nlhs, plhs );
    return;
  }
  
  
  /* Get input matrix dimensions */
  
  M = ( int )mxGetM( prhs[0] );
  N = ( int )mxGetN( prhs[0] );
  
  /* Check input dimensions: A vector of unity length is considered
   * invalid; matrices are invalid. */
  
  if( (M==1 && N==1) || (M!=1 && N!=1) ) {
    tfsaErr( "spec", "Input must be a vector" );
    winNTcheck( nlhs, plhs );
    return;
  }
  
  signal_length = (int)(M>N?M:N);
  
  
  /* Check second input */
  if( !GoodScalar( (MATRIX *)prhs[1] ) ) {
      tfsaErr( "spec", "Time resolution must be a scalar" );
      winNTcheck( nlhs, plhs );
      return;
  }
  else
    time_res = (int)*(mxGetPr( prhs[1]));
  
  /* Parameter value check */
  if( time_res < 1 ) {
    tfsaErr( "spec", "Time resolution must be greater than zero" );
    winNTcheck( nlhs, plhs );
    return;
  }
  
  if (time_res > signal_length) {
    tfsaErr( "spec", "Time resolution must be no greater than signal length" );
    winNTcheck( nlhs, plhs );
    return;
  }

  /* Check third input */
  if( !GoodScalar( (MATRIX *)prhs[2] ) ) {
    tfsaErr( "spec", "Smoothing window length must be a scalar" );
    winNTcheck( nlhs, plhs );
    return;
  }
  else
    window_length = (int)*(mxGetPr( prhs[2] ) );  /* Get value */
  
  
  /* Parameter value check */
  if( window_length < 1 ) {
    tfsaErr( "spec", "Window length must be greater than zero" );
    winNTcheck( nlhs, plhs );
    return;
  }
  
  /* Check window length against signal length */

  if( window_length > signal_length ) {
    tfsaWarning( "spec", "Window length has been truncated to signal length" );
    window_length = (int)signal_length;
  }
  
  /* Check for odd window length.. */
  if( fmod( window_length, 2) == 0 ){
    window_length = (int)(window_length - 1);
    tfsaWarning( "spec", "Window length is not odd, truncating to odd value." );
  } 

  
  /* Check fourth input*/
  if( MXSTRING( prhs[3] ) ) 
    {
      
      /* Get string from Matlab following p2-51 in Mex manual */
      
      n = (int)mxGetN( prhs[3] ) + 1;  /* string length + 1 for NULL */
      p = mxCalloc( n, sizeof( char ) );
      if( p == NULL ) 
	{
	  tfsaErr( "spec", "Internal memory allocation failure" );
	  winNTcheck( nlhs, plhs );
	  return;
	}
      
      if( mxGetString( prhs[3], p, n ) ) 
	{
	  tfsaErr( "spec", "Could not get smoothing window type" );
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
      else 
	{
	  tfsaErr( "spec", "Unknown smoothing window type" );
	  winNTcheck( nlhs, plhs );
	  return;
	}
    }
  else 
    {
      tfsaErr( "spec",  "Smoothing window type must be a string" );
      winNTcheck( nlhs, plhs );
      return;
    }

  /* Check fifth input if exist*/
  if( nrhs >= 5 ) {
    if( !GoodScalar( (MATRIX *)prhs[4] ) ) {
      tfsaErr( "spec", "FFT length must be a scalar" );
      winNTcheck( nlhs, plhs );
      return;
    }
    else
      fft_length = (int)*(mxGetPr( prhs[4]));
  }
  else
    fft_length = (int)signal_length;
  
  if( fft_length < 1 ) {
    tfsaErr( "spec", "FFT length must be greater than zero" );
    winNTcheck( nlhs, plhs );
    return;
  }
  
  if (fft_length < window_length) fft_length = window_length;


  if( nrhs == 6 ) {
    if( !GoodScalar( (MATRIX *)prhs[5] ) ) {
      tfsaErr( "spec", "STFT/Spectrogram flag must be a scalar" );
      winNTcheck( nlhs, plhs );
      return;
    }
    else
      stft_or_spec = abs( (int)*(mxGetPr( prhs[5])) );
  }
  else
    stft_or_spec = 0;
    

  
  signal_real = mxGetPr( prhs[0]);
  if (mxIsComplex(prhs[0]))
    signal_imag = mxGetPi( prhs[0]);
  else
    signal_imag = NULL;
  

  
  /* calculate radix-2 value equal or above fft_length */
  window_order = 0;
  window_r2 = 1;
  while( window_r2 < fft_length) {
      window_order++;
      window_r2 <<= 1;
  }


  /* Find out which is smaller, FFT length or signal length */
  N_min = (int) (window_r2 < signal_length ? window_r2 : signal_length );


  /* Calculate number of times for computing FFT. This
   * is equal to the number of estimated frequency points
   */
  numslices = (int) ceil((double)N_min/time_res);
  
  plhs[0] = MXCREATEFULL( (window_r2 / 2) + 1, numslices, COMPLEX );

  if( !plhs[0] ) {
      tfsaErr("spec", "Memory allocation failed" );
      winNTcheck( nlhs, plhs );
      return;
  }

  tfd_re = mxGetPr( plhs[0] );
  tfd_im = mxGetPi( plhs[0] );
  
  spec( signal_real, signal_imag, signal_length, tfd_re, tfd_im, numslices, 
        time_res, window_length, window_type, window_r2, window_order, stft_or_spec );
}


