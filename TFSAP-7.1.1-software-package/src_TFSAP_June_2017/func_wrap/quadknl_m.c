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
* gateway routine for quadknl.c
*************************************************/
#include <stdio.h>
#include <string.h>

#include "mex.h"
#include "matrix.h"

#include "arithm.h"		     /*	for complex vectors */
#include "tlocal.h"		     /* local function prototypes */
 
#include "quadknl.h"
#include "tfsa_c.h"

void mexFunction( int nlhs, MATRIX *plhs[], int nrhs, NEED_CONST MATRIX *prhs[] )
{
  int window_length;
  double *kernel_r;
  int full_kernel;
  double **G_real;
  int i, j, n;
  int kernel_type, window_type;
  double d_param, d_param_beta;
  int i_param;
  int size;
  char *p;

  /* Check number of input arguments */

  if( nrhs < 3 ) {
    tfsaErr( "quadknl", "Not enough input arguments" );
    winNTcheck( nlhs, plhs );
    return;
  }


  /* First input argument */

  if( MXSTRING( prhs[0] ) ) {

    /* Get string from Matlab following p2-51 in Mex manual */

    n = (int)mxGetN( prhs[0] ) + 1;  /* string length + 1 for NULL */
    p = mxCalloc( n, sizeof( char ) );
    if( p == NULL ) {
      tfsaErr( "quadknl", "Internal memory allocation failure" );
      winNTcheck( nlhs, plhs );
      return;
    }

    if( mxGetString( prhs[0], p, n ) ) {
      tfsaErr( "quadknl", "Could not get kernel_type" );
      winNTcheck( nlhs, plhs );
      return;
    }

    /* Check to see if valid string type:  */

    if( !strcmp( p, "wvd" ) || !strcmp( p, "WVD" ) )
      kernel_type = WVD;



    else if( !strcmp( p, "smoothed") || !strcmp( p, "SMOOTHED") )
      kernel_type = SMOOTHED;
    else if( !strcmp( p, "stft" ) || !strcmp( p, "STFT" ) )
      kernel_type = STFT;
    else if( !strcmp( p, "rm" ) || !strcmp( p, "RM" ) )
      kernel_type = RM;
    else if( !strcmp( p, "cw" ) || !strcmp( p, "CW" ) )
      kernel_type = CW;
    else if( !strcmp( p, "bjc") || !strcmp( p, "BJC") )
      kernel_type = BJC;
    else if( !strcmp( p, "zam") || !strcmp( p, "ZAM") )
      kernel_type = ZAM;
    else if( !strcmp( p, "b" ) || !strcmp( p, "B" ) )
      kernel_type = B;
    else if( !strcmp( p, "mb" ) || !strcmp( p, "MB" ) )
      kernel_type = MB;
    else if( !strcmp( p, "emb" ) || !strcmp( p, "EMB" ) )
      kernel_type = EMB;
    else {
      tfsaErr( "quadknl", "Unknown kernel type" );
      winNTcheck( nlhs, plhs );
      return;
    }
  }
  else {
    tfsaErr( "quadknl", "Kernel type must be a string" );
    winNTcheck( nlhs, plhs );
    return;
  }

  /* Check to see if too many input parameters have been given for the
   * distribution selected. */

  switch( kernel_type ) {

  case( WVD ):  /* These distributions require three input parameters */
  case( RM ):
  case( BJC ):

    if( nrhs > 3 ) {
      tfsaErr( "quadknl", "Too many input arguments" );
      winNTcheck( nlhs, plhs );
      return;
    }
    break;


  case( CW ):   /* These distributions require four input parameters */
  case( ZAM ):
  case( B ):
  case( MB ):

    if( nrhs > 4 ) {
      tfsaErr( "quadknl", "Too many input arguments" );
      winNTcheck( nlhs, plhs );
      return;
    }
    break;

  case( SMOOTHED ):  /* Smoothed WVD requires five inputs */
  case( STFT ):  

    if( nrhs > 5 ) {
      tfsaErr( "quadknl", "Too many input arguments" );
      winNTcheck( nlhs, plhs );
      return;
    }
    break;
  }


  /* Check second input parameter: window length  */

  if( !GoodScalar( (MATRIX *)prhs[1] ) ) {
    tfsaErr( "quadknl", "Window length must be a real scalar" );
    winNTcheck( nlhs, plhs );
    return;
  }

  window_length = (int)*( mxGetPr( prhs[1] ) );  /* get scalar value */

  if( window_length < 1 )  {
    tfsaErr( "quadknl", "Window length must be a non-negative scalar" );
    winNTcheck( nlhs, plhs );
    return;
  }

  if( !(window_length % 2) ) {
    tfsaErr( "quadknl", "Window length must be odd" );
    winNTcheck( nlhs, plhs );
    return;
  }


  /* Check third input parameter: full kernel flag  */

  if( !GoodScalar( (MATRIX *)prhs[2] ) ) {
    tfsaErr( "quadknl", "Full kernel flag must be either 0 or 1" );
    winNTcheck( nlhs, plhs );
    return;
  }

  full_kernel = (int)*mxGetPr( prhs[2] );      /* get value */

  if( (full_kernel != 0) && (full_kernel != 1) ) {
    tfsaErr( "quadknl", "Full kernel flag must be either 0 or 1" );
    winNTcheck( nlhs, plhs );
    return;
  }

  /* Validate input arguments */


  switch( kernel_type ) {

  case( SMOOTHED ):
  case( STFT ):

    if( nrhs < 5 ) {
      tfsaErr( "quadknl", "Smoothed kernel requires a smoothing window width and type" );
      winNTcheck( nlhs, plhs );
      return;
    }

    if( !GoodScalar( (MATRIX *)prhs[3] ) )  {
      tfsaErr( "quadknl", "Smoothing window length must be a scalar" );
      winNTcheck( nlhs, plhs );
      return;
    }

    i_param = (unsigned)*(mxGetPr( prhs[3] ));  /* window length */

    if( i_param < 1 )  {
      tfsaErr( "quadknl", "Smoothing window length must be non-negative" );
      winNTcheck( nlhs, plhs );
      return;
    }

    if( !(i_param % 2) ) {
      tfsaErr( "quadknl", "Smoothing window width must be odd" );
      winNTcheck( nlhs, plhs );
      return;
    }

    if( MXSTRING( prhs[4] ) ) {

      /* Get string from Matlab following p2-51 in Mex manual */

      n = (int)mxGetN( prhs[4] ) + 1;  /* string length + 1 for NULL */
      p = mxCalloc( n, sizeof( char ) );
      if( p == NULL ) {
	tfsaErr( "quadknl", "Internal memory allocation failure" );
	winNTcheck( nlhs, plhs );
	return;
      }

      if( mxGetString( prhs[4], p, n ) ) {
        tfsaErr( "quadknl", "Could not get smoothing window type" );
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
        tfsaErr( "quadknl", "Unknown smoothing window type" );
	winNTcheck( nlhs, plhs );
	return;
      }
    }
    else {
      tfsaErr( "quadknl", "Smoothing window type must be a string" );
      winNTcheck( nlhs, plhs );
      return;
    }
    break;

  case( RM ):

    if( full_kernel == 0 )  {
      full_kernel = 1;
      tfsaWarning( "quadknl", "Full kernel must be used with Rihaczek-Margenau distribution");
    }
    break;

  case( CW ):

    if( nrhs < 4 )  {
      tfsaErr( "quadknl", "Choi-Williams distribution requires a smoothing parameter (sigma)" );
      winNTcheck( nlhs, plhs );
      return;
    }


    if( !GoodScalar( (MATRIX *)prhs[3] ) )  {
      tfsaErr( "quadknl", "Smoothing parameter sigma must be a scalar" );
      winNTcheck( nlhs, plhs );
      return;
    }

    d_param = *(mxGetPr( prhs[3] ));       /* Get sigma value */

    if( d_param < 0.0 )  {
      tfsaErr( "quadknl", "Smoothing parameter sigma must be non-negative" );
      winNTcheck( nlhs, plhs );
      return;
    }

    break;

  case( ZAM ):
    if( nrhs < 4 ) {
      tfsaErr( "quadknl", "ZAM distribution requires the parameter 'a'" );
      winNTcheck( nlhs, plhs );
      return;
    }


    if( !GoodScalar( (MATRIX *)prhs[3] ) )  {
      tfsaErr( "quadknl", "ZAM distribution parameter 'a' must be a scalar" );
      winNTcheck( nlhs, plhs );
      return;
    }

    d_param = *(mxGetPr( prhs[3]));

    if( d_param < 0.0 )  {
      tfsaErr( "quadknl", "ZAM distribution parameter 'a' must be non-negative" );
      winNTcheck( nlhs, plhs );
      return;
    }

    if( d_param < 1 ) {
      tfsaErr( "quadknl", "ZAM distribution parameter 'a' must be equal to or greater than one" );
      winNTcheck( nlhs, plhs );
      return;
    }
    break;

  case( B ):

    if( nrhs < 4 )  {
      tfsaErr( "quadknl", "B  distribution requires a smoothing parameter (beta)" );
      winNTcheck( nlhs, plhs );
      return;
    }


    if( !GoodScalar( (MATRIX *)prhs[3] ) )  {
      tfsaErr( "quadknl", "Smoothing parameter beta must be a scalar" );
      winNTcheck( nlhs, plhs );
      return;
    }

    d_param = *(mxGetPr( prhs[3] ));       /* Get beta value */

    if( d_param < 0.0 || d_param > 1.0)  {
      tfsaErr( "quadknl", "Smoothing parameter beta must be between 0 and 1" );
      winNTcheck( nlhs, plhs );
      return;
    }
    break;  

  case( MB ):

    if( nrhs < 4 )  {
      tfsaErr( "quadknl", "Modified B-distribution requires a smoothing parameter (alpha)" );
      winNTcheck( nlhs, plhs );
      return;
    }


    if( !GoodScalar( (MATRIX *)prhs[3] ) )  {
      tfsaErr( "quadknl", "Smoothing parameter alpha must be a scalar" );
      winNTcheck( nlhs, plhs );
      return;
    }

    d_param = *(mxGetPr( prhs[3] ));       /* Get alpha value */

    if( d_param < 0.0 || d_param > 1.0)  {
      tfsaErr( "quadknl", "Smoothing parameter alpha must be between 0 and 1" );
      winNTcheck( nlhs, plhs );
      return;
    }

    break;


    case( EMB ):

    if( nrhs < 5 )  {
      tfsaErr( "quadknl", "Extended Modified B-distribution requires two smoothing parameters (alpha, beta)" );
      winNTcheck( nlhs, plhs );
      return;
    }


    if( !GoodScalar( (MATRIX *)prhs[3] ) )  {
      tfsaErr( "quadknl", "Smoothing parameter alpha must be a scalar" );
      winNTcheck( nlhs, plhs );
      return;
    }

    d_param = *(mxGetPr( prhs[3] ));       /* Get alpha value */

    if( d_param < 0.0 || d_param > 1.0)  {
      tfsaErr( "quadknl", "Smoothing parameter alpha must be between 0 and 1" );
      winNTcheck( nlhs, plhs );
      return;
    }

    d_param_beta = *(mxGetPr( prhs[4] ));       /* Get alpha beta */
    
    if( d_param_beta < 0.0 || d_param_beta > 1.0)  {
        tfsaErr( "quadknl", "Smoothing parameter beta must be between 0 and 1" );
        winNTcheck( nlhs, plhs );
        return;
    }

    break;
  }  /* switch( kernel_type ) */


  /* kernel options */

  if( full_kernel )
    size = window_length;
  else
    size = (window_length+1)/2;


  /* Check length of smoothed window for smoothed Wigner-Ville */

  if(( kernel_type == SMOOTHED ) || ( kernel_type == STFT )){
    if( i_param > window_length ) {
      i_param = window_length;
      tfsaWarning("quadknl","Smoothing window width has been truncated to window length");
    }
  }



  /* check output arguments */

  if( nlhs > 1 ) {
    tfsaErr( "quadknl", "Too many output arguments" );
    winNTcheck( nlhs, plhs );
    return;
  }


  /* allocate output matrix, and check memory allocation */

  plhs[0] = MXCREATEFULL( size, size, REAL);
  if( plhs[0] == NULL ) {
    tfsaErr( "quadknl", "Internal memory allocation failure: kernel_r" );
    winNTcheck( nlhs, plhs );
    return;
  }

  kernel_r = mxGetPr( plhs[0] );

  G_real = (double **)mxCalloc( size, sizeof(double *) );

  if( G_real == NULL ) {
    tfsaErr( "quadknl", "Internal memory allocation failure: G_real" );
    winNTcheck( nlhs, plhs );
    return;
  }

  G_real[0] = (double *) mxCalloc( (int)(size*size), sizeof(double));

  if( G_real[0] == NULL) {
    tfsaErr( "quadknl", "Internal memory allocation failure: G_real[0]" );
    winNTcheck( nlhs, plhs );
    return;
  }

  /* Set up pointers for pointer array */
  for(i = 1; i < size; i++ )
    G_real[i] = G_real[0] + (unsigned long)size * i;


  /* Call respective distributions */

  switch( kernel_type ) {
  case( WVD ):
    wvd_kernel( G_real, window_length, full_kernel );
    break;

  case( SMOOTHED ):
    smoothedwvd_kernel( G_real, window_length, full_kernel, i_param, window_type );
    break;

  case( STFT ):
    stft_kernel( G_real, window_length, full_kernel, i_param, window_type );
    break;

  case( RM ):
    complex_rm_kernel( G_real, window_length );
    break;

  case( CW ):
    cw_kernel( G_real, window_length, full_kernel, d_param );
    break;

  case( BJC ):
    bjc_kernel( G_real, window_length, full_kernel );
    break;

  case( ZAM ):
    zam_kernel( G_real, window_length, full_kernel, d_param );
    break;

  case( B ):
    b_kernel( G_real, window_length, full_kernel, d_param );
    break;

  case( MB ):
    mb_kernel( G_real, window_length, full_kernel, d_param );
    break;
    
  case EMB:
      emb_kernel( G_real, window_length, full_kernel, d_param, d_param_beta);
	  break;
  }

  for( i = 0; i < size; i++ ) {
    for( j = 0; j< size; j++ ) {
      kernel_r[i + j * (unsigned long)size] = G_real[i][j];
    }
  }
}










