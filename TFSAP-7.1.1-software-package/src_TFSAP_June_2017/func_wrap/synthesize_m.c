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
* gateway routine for synthesize.c.
*************************************************/

#include <math.h>
#include "mex.h"
#include <matrix.h>
#include <string.h>

#include "arithm.h"		     /* routines for complex data */
#include "tlocal.h"		     /* local function prototypes */

#include "synthesize.h"
#include "tfsa_c.h"
#include "window.h"



void mexFunction( int nlhs, MATRIX *plhs[], int nrhs, NEED_CONST MATRIX *prhs[] )
{
  double *signal_real = NULL, *signal_imag = NULL;
  double *synth_signal_re, *synth_signal_im;
  int window_length;
  int window_type;
  int signal_length;
  int M,N,M_y,N_y;
  int n;
  int analysis_type;
  double tol;
  char *p;
  double *tfd_re = NULL, *tfd_im = NULL;
  double *rankone_err = NULL;


  /* TMP !!!! */
  /* char db_msg[50]; */




  /* basic input--output number check */
  
  if( nrhs < 3 ) {
    tfsaErr( "synthesize", "Not enough input arguments" );
    winNTcheck( nlhs, plhs );
    return;
  }
  
  if( nrhs > 6 ) {
    tfsaErr( "synthesize", "Too many input arguments" );
    winNTcheck( nlhs, plhs );
    return;
  }
  
  if( nlhs > 2 ) {
    tfsaErr( "synthesize", "Too many output arguments" );
    winNTcheck( nlhs, plhs );
    return;
  }

  
  
  /* Check first input */
  
  if( !mxIsNumeric( prhs[0] ) || !MXFULL( prhs[0] ) ) {
    tfsaErr( "synthesize", "Input must be a vector" );
    winNTcheck( nlhs, plhs );
    return;
  }
  
  
  /* Get input matrix dimensions */
  M_y = ( int )mxGetM( prhs[0] );
  N_y = ( int )mxGetN( prhs[0] );
  
  /* Check input dimensions: A vector of unity length is considered invalid. */
  if( N_y == 1 )  {
    tfsaErr( "synthesize", "Input TFD must be a matrix" );
    winNTcheck( nlhs, plhs );
    return;
  }

  tfd_re = mxGetPr( prhs[0] );
  tfd_im = mxGetPi( prhs[0] );


  signal_length = (int)N_y;  

  /* Check second input*/
  if( MXSTRING( prhs[1] ) ) 
    {
      
      n = (int)mxGetN( prhs[1] ) + 1;  /* string length + 1 for NULL */
      p = mxCalloc( n, sizeof( char ) );
      if( p == NULL ) 
	{
	  tfsaErr( "synthesize", "Internal memory allocation failure" );
	  winNTcheck( nlhs, plhs );
	  return;
	}
      
      if( mxGetString( prhs[1], p, n ) ) 
	{
	  tfsaErr( "synthesize", "Could not get analysis type" );
	  winNTcheck( nlhs, plhs );
	  return;
	}
      
      if( !strcmp( p, "idft" ) || !strcmp( p, "IDFT" ) )
	analysis_type = IDFT;
      else if( !strcmp( p, "ola" ) || !strcmp( p, "OLA" ) )
	analysis_type = OLA;
      else if( !strcmp( p, "mstft" ) || !strcmp( p, "MSTFT" ) )
	analysis_type = MSTFT;
      else if( !strcmp( p, "mspec" ) || !strcmp( p, "MSPEC" ) )
	analysis_type = MSPEC;
      else if( !strcmp( p, "wvd" ) || !strcmp( p, "WVD" ) )
	analysis_type = MWVD;
      else 
	{
	  tfsaErr( "synthesize", "Unknown analysis type" );
	  winNTcheck( nlhs, plhs );
	  return;
	}

    } else {

      tfsaErr( "synthesize",  "Analysis type must be a string" );
      winNTcheck( nlhs, plhs );
      return;
    }

  
  /* Check third input */
  if( !GoodScalar( (MATRIX *)prhs[2] ) ) {
    tfsaErr( "synthesize", "Smoothing window length must be a scalar" );
    winNTcheck( nlhs, plhs );
    return;
  }
  else
    window_length = (int)*(mxGetPr( prhs[2] ) );  /* Get value */
  
  
  /* Parameter value check */
  if( window_length < 1 ) {
    tfsaErr( "synthesize", "Window length must be greater than zero" );
    winNTcheck( nlhs, plhs );
    return;
  }
  
  /* Check window length against signal length */

  if( window_length > signal_length ) {
    tfsaWarning( "synthesize", "Window length has been truncated to signal length" );
    window_length = (int)signal_length;
  }
  
  /* Check for odd window length.. */
  if( fmod( window_length, 2) == 0 ){
    window_length = (int)(window_length - 1);
    tfsaWarning( "synthesize", "Window length is not odd, truncating to odd value." );
  } 
  
    
  /* Check fourth input*/

  if( analysis_type == MWVD ){
    /* Original signal, used to correct phase */
    if( nrhs == 4 ){
      if( !mxIsNumeric( prhs[3] ) || !MXFULL( prhs[3] ) ) {
	tfsaErr( "synthesize", "Input must be a vector" );
	winNTcheck( nlhs, plhs );
	return;
      }
  
      /* Get input matrix dimensions */
      M = ( int )mxGetM( prhs[3] );
      N = ( int )mxGetN( prhs[3] );
  
      /* Check input dimensions: A vector of unity length is considered
       * invalid; matrices are invalid. */
  
      if( (M==1 && N==1) || (M!=1 && N!=1) ) {
	tfsaErr( "synthesize", "Input must be a vector" );
	winNTcheck( nlhs, plhs );
	return;
      }

      signal_real = mxGetPr( prhs[3]);
      if( mxIsComplex(prhs[3]) )
	signal_imag = mxGetPi( prhs[3] );
    }

  } else {


    if( MXSTRING( prhs[3] ) ) 
    {
      
      /* Get string from Matlab following p2-51 in Mex manual */
      
      n = (int)mxGetN( prhs[3] ) + 1;  /* string length + 1 for NULL */
      p = mxCalloc( n, sizeof( char ) );
      if( p == NULL ) 
	{
	  tfsaErr( "synthesize", "Internal memory allocation failure" );
	  winNTcheck( nlhs, plhs );
	  return;
	}
      
      if( mxGetString( prhs[3], p, n ) ) 
	{
	  tfsaErr( "synthesize", "Could not get smoothing window type" );
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
	  tfsaErr( "synthesize", "Unknown smoothing window type" );
	  winNTcheck( nlhs, plhs );
	  return;
	}

    } else {
      tfsaErr( "synthesize",  "Smoothing window type must be a string" );
      winNTcheck( nlhs, plhs );
      return;
    }
  }
    
  /* Check fifth input if exist*/
  if( nrhs > 4 ){
    if( !GoodScalar( (MATRIX *)prhs[4] ) ) {
      tfsaErr( "synthesize", "Tolerance value must be a scalar" );
      winNTcheck( nlhs, plhs );
      return;
    }
    else
      tol = (double)*(mxGetPr( prhs[4]));
  } else {
    tol = TOL_DEFAULT;
  }


  /* Check sixth input if exist*/
  if( nrhs == 6 ){
    /* Original signal, used to correct phase */

      if( !mxIsNumeric( prhs[5] ) || !MXFULL( prhs[5] ) ) {
	tfsaErr( "synthesize", "Input must be a vector" );
	winNTcheck( nlhs, plhs );
	return;
      }
  
      /* Get input matrix dimensions */
      M = ( int )mxGetM( prhs[5] );
      N = ( int )mxGetN( prhs[5] );
  
      /* Check input dimensions: A vector of unity length is considered
       * invalid; matrices are invalid. */
  
      if( (M==1 && N==1) || (M!=1 && N!=1) ) {
	tfsaErr( "synthesize", "Input must be a vector" );
	winNTcheck( nlhs, plhs );
	return;
      }

      signal_real = mxGetPr( prhs[5]);
      if( mxIsComplex(prhs[5]) )
	signal_imag = mxGetPi( prhs[5] );
  }



  /* Allocate memory and space for the returned synthesized signal */
  plhs[0] = MXCREATEFULL (signal_length, 1, COMPLEX);
  
  if (plhs[0] == NULL){
      tfsaErr("synthesize", "Memory allocation failed" );
      winNTcheck( nlhs, plhs );
      return;
  }
  synth_signal_re = mxGetPr( plhs[0] );
  synth_signal_im = mxGetPi( plhs[0] );

  
  if( analysis_type == MWVD ){
    plhs[1] = MXCREATEFULL (1, 1, REAL);
  
    if (plhs[1] == NULL){
      tfsaErr("synthesize", "Memory allocation failed" );
      winNTcheck( nlhs, plhs );
      return;
    }
    rankone_err = mxGetPr( plhs[1] );
  }



  synthesize( analysis_type, synth_signal_re, synth_signal_im, signal_length, 
	      window_length, window_type, tol, signal_real, signal_imag, 
	      tfd_re, tfd_im, M_y, N_y, rankone_err );

}


