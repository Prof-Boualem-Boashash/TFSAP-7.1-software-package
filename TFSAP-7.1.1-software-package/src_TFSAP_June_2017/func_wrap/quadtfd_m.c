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
* gateway routine for quadtfd.c
*************************************************/

#include <string.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"

#include "arithm.h"		     /* routines for complex data */
#include "tlocal.h"		     /* local function prototypes */

#include "quadtfd.h"
#include "quadknl.h"
#include "tfsa_c.h"

void mexFunction( int nlhs, MATRIX *plhs[], int nrhs, NEED_CONST MATRIX *prhs[])
{
  double *result_r, *result_i;
  double *signal_r, *signal_i;
  int window_length, time_res;  // changes i made
  int signal_length;
  int M,N;
  int num_slices, window_order, window_r2;
  int i, j;
  double **G_real, **G_imag;
  int kernel_type, window_type;
  int size;
  int fft_length;

  int size_to_check;
  
  char *p;
  int m,n;

  int i_param;
  double d_param;
  double d_param_beta;

  
  /* Basic input--output argument number check */

  if( nrhs < 4 ) {
    tfsaErr( "quadtfd", "Not enough input arguments" );
    winNTcheck( nlhs, plhs );
    return;
  }

  if( nlhs > 1 ) {
    tfsaErr( "quadtfd", "Too many output arguments" );
    winNTcheck( nlhs, plhs );
    return;
  }


  /* Check and parse input variables */

  /* First input -- can be complex */

  if( !mxIsNumeric( prhs[0] ) || !MXFULL( prhs[0] ) ) {
    tfsaErr( "quadtfd", "Input signal must be a vector" );
    winNTcheck( nlhs, plhs );
    return;
  }

  M = (int)mxGetM( prhs[0] );
  N = (int)mxGetN( prhs[0] );

  if( (M==1 && N==1) || (M!=1 && N!=1) ) {
    tfsaErr( "quadtfd", "Input signal must be a vector" );
    winNTcheck( nlhs, plhs );
    return;
  }

  /* Choose the larger row or column dimension for signal length */

  signal_length = (unsigned)( M>N ? M:N );

  /* Dereference input */

  signal_r = mxGetPr( prhs[0] );
  if( mxIsComplex( prhs[0] ) )
    signal_i = mxGetPi( prhs[0]);
  else
    signal_i = NULL;

  /* Check second input parameter: window length  */

  if( !GoodScalar( (MATRIX *)prhs[1] ) ) {
    tfsaErr( "quadtfd", "Window length must be a real scalar" );
    winNTcheck( nlhs, plhs );
    return;
  }

  window_length = (int)*( mxGetPr( prhs[1] ) );  /* get scalar value */

  if( window_length < 1 )  {
    tfsaErr( "quadtfd", "Window length must be greater than zero" );
    winNTcheck( nlhs, plhs );
    return;
  }

  if( window_length > signal_length  ) {
    tfsaWarning( "quadtfd", "Window length has been truncated to signal length" );
    window_length = signal_length;

    if( !(window_length%2) )   /* Make window length odd if even */
      window_length--;
  }

   

  if( !(window_length % 2) ) {
    tfsaErr( "quadtfd", "Window length must be odd" );
    winNTcheck( nlhs, plhs );
    return;
  }

  /* Check third input parameter: time res */

  if( !GoodScalar( (MATRIX *)prhs[2] ) ) {
    tfsaErr( "quadtfd", "Time resolution must be a real scalar" );
    winNTcheck( nlhs, plhs );
    return;
  }

  time_res = (int)*(mxGetPr( prhs[2]));

  if( time_res < 1 )  {
    tfsaErr( "quadtfd", "Time resolution must be greater than zero" );
    winNTcheck( nlhs, plhs );
    return;
  }

  if( time_res > signal_length )  {
    tfsaWarning( "quadtfd", "Time resolution set to signal length" );
    time_res = (int)signal_length;

  }



  /* calculate distribution size */

  num_slices = (int) ceil((double)signal_length/time_res);


  /* Fourth input argument */

  if( MXSTRING( prhs[3] ) ) {

    /* Get string from Matlab following p2-51 in Mex manual */

    n = (int)mxGetN( prhs[3] ) + 1;  /* string length + 1 for NULL */
    p = mxCalloc( n, sizeof( char ) );
    if( p == NULL ) {
      tfsaErr( "quadtfd", "Internal memory allocation failure" );
      winNTcheck( nlhs, plhs );
      return;
    }

    if( mxGetString( prhs[3], p, n ) ) {
      tfsaErr( "quadtfd", "Could not get kernel type" );
      winNTcheck( nlhs, plhs );
      return;
    }

    if( !strcmp( p, "wvd" ) || !strcmp( p, "WVD" ) )
      kernel_type = WVD;

    else if( !strcmp( p, "smoothed" ) || !strcmp( p, "SMOOTHED" ) )
      kernel_type = SMOOTHED;
    else if( !strcmp( p, "specx" ) ||  !strcmp( p, "SPECX" ) )
      kernel_type = STFT;
    else if( !strcmp( p, "rm" ) || !strcmp( p, "RM" ))
      kernel_type = RM;
    else if( !strcmp( p, "cw" ) || !strcmp( p, "CW" ) )
      kernel_type = CW;
    else if( !strcmp( p, "bjc" ) || !strcmp( p, "BJC" ) )
      kernel_type = BJC;
    else if( !strcmp( p, "zam" ) || !strcmp( p, "ZAM" ) )
      kernel_type = ZAM;
    else if( !strcmp( p, "b" ) || !strcmp( p, "B" ) )
      kernel_type = B;
    else if( !strcmp( p, "mb" ) || !strcmp( p, "MB" ) )
      kernel_type = MB;
	else if( !strcmp( p, "emb" ) || !strcmp( p, "EMB" ) )
      kernel_type = EMB;
    else {
      tfsaErr( "quadtfd", "Unknown kernel type" );
      winNTcheck( nlhs, plhs );
      return;
    }
  }
  else {
    tfsaErr( "quadtfd", "Smoothing window type must be a string" );
    winNTcheck( nlhs, plhs );
    return;
    }



  if( (kernel_type != USER) && (kernel_type != RM) ) {

    /* using a real, symmetric kernel */

    size = (window_length+1)/2;

    /* check additional options */

    switch( kernel_type ) {

    case( SMOOTHED ):  /* Smoothed WVD requires six inputs, with optional fft_length */
    case( STFT ):      /* Now STFT             */

      /* check number of input parameters */

      if( nrhs < 6 ) {
	tfsaErr( "quadtfd", "Not enough input arguments" );
	winNTcheck( nlhs, plhs );
	return;
      }

      if( nrhs > 7 ) {
	tfsaErr( "quadtfd", "Too many input arguments" );
	winNTcheck( nlhs, plhs );
	return;
      }


      /* Process fifth input */

      if( !GoodScalar( (MATRIX *)prhs[4] ) )  {
	tfsaErr( "quadtfd", "Smoothing window length must be a scalar" );
	winNTcheck( nlhs, plhs );
	return;
      }

      i_param = (int)*(mxGetPr( prhs[4] ));

      if( i_param < 1 )  {
	tfsaErr( "quadtfd", "Smoothing window length must be greater than zero" );
	winNTcheck( nlhs, plhs );
	return;
      }

      if( i_param > window_length  )  {
	tfsaWarning( "quadtfd", "Smoothing window length has been truncated to window length");
	i_param = window_length;

	if( !(i_param%2) )   /* Make smoothing window length odd if even */
	  i_param--;
      }

      if( !(i_param % 2) ) {
	tfsaErr( "quadtfd", "Smoothing window width must be odd" );
	winNTcheck( nlhs, plhs );
	return;
      }

      /* Process sixth input */

      if( MXSTRING( prhs[5] ) ) {

	/* Get string from Matlab following p2-51 in Mex manual */

	n = (int)mxGetN( prhs[5] ) + 1;  /* string length + 1 for NULL */
	p = mxCalloc( n, sizeof( char ) );
	if( p == NULL ) {
	  tfsaErr( "quadtfd", "Internal memory allocation failure" );
	  winNTcheck( nlhs, plhs );
	  return;
	}

	if( mxGetString( prhs[5], p, n ) ) {
	  tfsaErr( "quadtfd", "Could not get smoothing window type" );
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
	  tfsaErr( "quadtfd", "Unknown smoothing window type" );
	  winNTcheck( nlhs, plhs );
	  return;
	}
      }
      else {
	tfsaErr( "quadtfd",  "Smoothing window type must be a string" );
	winNTcheck( nlhs, plhs );
	return;
      }

      /* Process seventh parameter if given -- fft length */

      if( nrhs == 7 ) {
	if( !GoodScalar( (MATRIX *)prhs[6] ) )  {
	  tfsaErr( "quadtfd", "FFT length must be a real scalar" );
	  winNTcheck( nlhs, plhs );
	  return;
	}

	fft_length = (int)*(mxGetPr( prhs[6] ));   /* get value */

	if( fft_length < 0 )  {
	  tfsaErr( "quadtfd", "FFT length must be greater than zero" );
	  winNTcheck( nlhs, plhs );
	  return;
	}
      }
      else
	fft_length = window_length; 	/* default */

      break;

    case( CW ):

      if( nrhs < 5 )  {
	tfsaErr( "quadtfd", "Not enough input arguments" );
	winNTcheck( nlhs, plhs );
	return;
      }

      if( nrhs > 6 ) {
	tfsaErr( "quadtfd", "Too many input arguments" );
	winNTcheck( nlhs, plhs );
      return;
      }


      /* Process fifth input parameter */

      if( !GoodScalar( (MATRIX *)prhs[4] ) )  {
	tfsaErr( "quadtfd", "Smoothing parameter sigma must be a scalar" );
	winNTcheck( nlhs, plhs );
	return;
      }

      d_param = *(mxGetPr( prhs[4] ));		/* Get sigma value */

      if( d_param < 0.0 )  {
	tfsaErr( "quadknl", "Smoothing parameter sigma must be non-negative" );
	winNTcheck( nlhs, plhs );
	return;
      }

      /* Process sixth input parameter if given */

      if( nrhs == 6 ) {			       /* fft length option ? */
	if( !GoodScalar( (MATRIX *)prhs[5] ) )  {
	  tfsaErr( "quadtfd", "FFT length must be a real scalar" );
	  winNTcheck( nlhs, plhs );
	  return;
	}

	fft_length = (int)*(mxGetPr( prhs[5] ));   /* get value */

	if( fft_length < 0 )  {
	  tfsaErr( "quadtfd", "FFT length must be greater than zero" );
	  winNTcheck( nlhs, plhs );
	  return;
	}
      }
      else
	fft_length = window_length;		/* default */

      break;

    case( ZAM ):

      if( nrhs < 5 )  {
	tfsaErr( "quadtfd", "Not enough input arguments" );
	winNTcheck( nlhs, plhs );
	return;
      }

      if( nrhs > 6 ) {
	tfsaErr( "quadtfd", "Too many input arguments" );
	winNTcheck( nlhs, plhs );
      return;
      }


      /* Process fifth input parameter */

      if( !GoodScalar( (MATRIX *)prhs[4] ) )  {
	tfsaErr( "quadtfd", "ZAM parameter 'a' must be a scalar" );
	winNTcheck( nlhs, plhs );
	return;
      }

      d_param = *(mxGetPr( prhs[4] ));		/* Get sigma value */

      if( d_param < 0.0	)  {
	tfsaErr( "quadknl", "ZAM parameter 'a' must be non-negative" );
	winNTcheck( nlhs, plhs );
	return;
      }

      /* Process sixth input parameter if exist */

      if( nrhs == 6 ) {			       /* fft length option ? */
	if( !GoodScalar( (MATRIX *)prhs[5] ) )  {
	  tfsaErr( "quadtfd", "FFT length must be a real scalar" );
	  winNTcheck( nlhs, plhs );
	  return;
	}

	fft_length = (int)*(mxGetPr( prhs[5] ));   /* get value	*/

	if( fft_length < 0 )  {
	  tfsaErr( "quadtfd", "FFT length must be greater than zero" );
	  winNTcheck( nlhs, plhs );
	  return;
	}
      }
      else
	fft_length = window_length;		/* default */

      break;

      case( B ):

      if( nrhs < 5 )  {
	tfsaErr( "quadtfd", "Not enough input arguments" );
	winNTcheck( nlhs, plhs );
	return;
      }

      if( nrhs > 6 ) {
	tfsaErr( "quadtfd", "Too many input arguments" );
	winNTcheck( nlhs, plhs );
      return;
      }

      /* Process fifth input parameter */

      if( !GoodScalar( (MATRIX *)prhs[4] ) )  {
	tfsaErr( "quadtfd", " Parameter beta must be a scalar" );
	winNTcheck( nlhs, plhs );
	return;
      }

      d_param = *(mxGetPr( prhs[4] ));		/* Get beta value */

      if( d_param < 0.0 || d_param > 1.0)  {
	tfsaErr( "quadknl", "Parameter beta must lie between 0 and 1." );
	winNTcheck( nlhs, plhs );
	return;
      }

      /* Process sixth input parameter if given */

      if( nrhs == 6 ) {			       /* fft length option ? */
	if( !GoodScalar( (MATRIX *)prhs[5] ) )  {
	  tfsaErr( "quadtfd", "FFT length must be a real scalar" );
	  winNTcheck( nlhs, plhs );
	  return;
	}

	fft_length = (int)*(mxGetPr( prhs[5] ));   /* get value */

	if( fft_length < 0 )  {
	  tfsaErr( "quadtfd", "FFT length must be greater than zero" );
	  winNTcheck( nlhs, plhs );
	  return;
	}
      }
      else
	fft_length = window_length;		/* default */

      break;

      case( MB ):

      if( nrhs < 5 )  {
	tfsaErr( "quadtfd", "Not enough input arguments" );
	winNTcheck( nlhs, plhs );
	return;
      }

      if( nrhs > 6 ) {
	tfsaErr( "quadtfd", "Too many input arguments" );
	winNTcheck( nlhs, plhs );
      return;
      }

      /* Process fifth input parameter */

      if( !GoodScalar( (MATRIX *)prhs[4] ) )  {
	tfsaErr( "quadtfd", " Parameter alpha must be a scalar" );
	winNTcheck( nlhs, plhs );
	return;
      }

      d_param = *(mxGetPr( prhs[4] ));		/* Get alpha value */

      if( d_param < 0.0 || d_param > 1.0)  {
	tfsaErr( "quadknl", "Parameter alpha must lie between 0 and 1" );
	winNTcheck( nlhs, plhs );
	return;
      }

      /* Process sixth input parameter if given */

      if( nrhs == 6 ) {			       /* fft length option ? */
	if( !GoodScalar( (MATRIX *)prhs[5] ) )  {
	  tfsaErr( "quadtfd", "FFT length must be a real scalar" );
	  winNTcheck( nlhs, plhs );
	  return;
	}

	fft_length = (int)*(mxGetPr( prhs[5] ));   /* get value */

	if( fft_length < 0 )  {
	  tfsaErr( "quadtfd", "FFT length must be greater than zero" );
	  winNTcheck( nlhs, plhs );
	  return;
	}
      }
      else
	fft_length = window_length;		/* default */

      break;




//  EMBD Parameter Settings -----------------


	  case( EMB ):

      if( nrhs < 6 )  {
	tfsaErr( "quadtfd", "Not enough input arguments" );
	winNTcheck( nlhs, plhs );
	return;
      }

      if( nrhs > 7 ) {
	tfsaErr( "quadtfd", "Too many input arguments" );
	winNTcheck( nlhs, plhs );
      return;
      }

      /* Process fifth input parameter */

      if( !GoodScalar( (MATRIX *)prhs[4] ) )  {
	tfsaErr( "quadtfd", " Parameter alpha must be a scalar" );
	winNTcheck( nlhs, plhs );
	return;
      }

	  if( !GoodScalar( (MATRIX *)prhs[5] ) )  {
	tfsaErr( "quadtfd", " Parameter beta must be a scalar" );
	winNTcheck( nlhs, plhs );
	return;
      }

      d_param = *(mxGetPr( prhs[4] ));		/* Get alpha value */

	  d_param_beta = *(mxGetPr( prhs[5] ));		/* Get beta value */

      if( d_param < 0.0 || d_param > 1.0)  {
	tfsaErr( "quadknl", "Parameter alpha must lie between 0 and 1" );
	winNTcheck( nlhs, plhs );
	return;
      }

	  if( d_param_beta < 0.0 || d_param_beta > 1.0)  {
	tfsaErr( "quadknl", "Parameter alpha must lie between 0 and 1" );
	winNTcheck( nlhs, plhs );
	return;
      }

      /* Process sixth input parameter if given */

      if( nrhs == 7 ) {			       /* fft length option ? */
	if( !GoodScalar( (MATRIX *)prhs[6] ) )  {
	  tfsaErr( "quadtfd", "FFT length must be a real scalar" );
	  winNTcheck( nlhs, plhs );
	  return;
	}

	fft_length = (int)*(mxGetPr( prhs[6] ));   /* get value */

	if( fft_length < 0 )  {
	  tfsaErr( "quadtfd", "FFT length must be greater than zero" );
	  winNTcheck( nlhs, plhs );
	  return;
	}
      }
      else
	fft_length = window_length;		/* default */

      break;

// END EMBD Parameter Settings -----------------

    case( WVD ):   /* requires four input parameters */
    case( RM ):
    case( BJC ):

      if( nrhs > 5 ) {
	tfsaErr( "quadtfd", "Too many input arguments" );
	winNTcheck( nlhs, plhs );
	return;
      }

      /* Process fifth input parameter if given */

      if( nrhs == 5 ) {			       /* fft length option ? */
	if( !GoodScalar( (MATRIX *)prhs[4] ) )  {
	  tfsaErr( "quadtfd", "FFT length must be a real scalar" );
	  winNTcheck( nlhs, plhs );
	  return;
	}

	fft_length = (int)*(mxGetPr( prhs[4] ));   /* get value */

	if( fft_length < 0 )  {
	  tfsaErr( "quadtfd", "FFT length must be greater than zero" );
	  winNTcheck( nlhs, plhs );
	  return;
	}
      }
      else
	  fft_length = window_length;		  /* default */

      break;

    }  /* end of switch( kernel_type ) */


    if (fft_length < window_length) fft_length = window_length;

    /* calculate radix-2 value equal or above window_length */
    window_order = 0;
    window_r2 = 1;
    while( window_r2 < fft_length) {
      window_order++;
      window_r2 <<= 1;
    }

    plhs[0] = MXCREATEFULL( window_r2, num_slices, REAL );
    if (plhs[0] == NULL) {
      tfsaErr( "quadtfd", "Memory allocation failure: plhs");
      winNTcheck( nlhs, plhs );
      return;
    }
    result_r = mxGetPr( plhs[0]);

    G_real = (double **)mxCalloc( size, sizeof(double *));
    if (G_real==NULL) {
      tfsaErr( "quadtfd", "Memory allocation failed: G_real");
      winNTcheck( nlhs, plhs );
      return;
    }
    G_real[0] = (double *)mxCalloc( (int)(size*size), sizeof(double));
    if (G_real[0]==NULL) {
      tfsaErr( "quadtfd", "Memory allocation failure: G_real[0]");
      winNTcheck( nlhs, plhs );
      return;
    }
    for (i = 1; i < size; i++)
      G_real[i] = G_real[0] + (unsigned long)size * i;

    /* generate kernel */
	/*mexPrintf("\n\nThe kernel is = %d \n", kernel_type);*/

    switch( kernel_type) {
    case WVD:
	  wvd_kernel( G_real, window_length, 0);
     
      break;

    case SMOOTHED:
      
      smoothedwvd_kernel( G_real, window_length, 0, i_param, window_type);
	 
      break;
    case STFT:
      stft_kernel( G_real, window_length, 0, i_param, window_type);
     
	  break;
    case CW:
      cw_kernel( G_real, window_length, 0, d_param);
     
	  break;
    case BJC:
      bjc_kernel( G_real, window_length, 0);
     
	  break;
    case ZAM:
      zam_kernel( G_real, window_length, 0, d_param);
     
	  break;
    case B:
	  b_kernel( G_real, window_length, 0, d_param);
	 
	  break;
    case MB:
      mb_kernel( G_real, window_length, 0, d_param);
	 
	  break;
    case EMB:
      emb_kernel( G_real, window_length, 0, d_param, d_param_beta);
	 
	  break;
    }


//mexPrintf("\n\nThe value of size is = %d", size);
//mexPrintf("\n\nThe value of window_length is = %u \n", analytic_or_real);



    quad_realG_symmG( signal_r, signal_i, signal_length, result_r,
		      num_slices, time_res, window_length,
		      window_r2, window_order, G_real);   // changes i made 

    mxFree( (void *)(G_real[0]));
    mxFree( (void *)G_real);
  }

  else if (kernel_type == RM) {

    size = window_length;

    if (nrhs < 5)
      fft_length = window_length;
    else{
      fft_length = (unsigned)*(mxGetPr( prhs[4]));
    }

    if (fft_length < window_length) fft_length = window_length;

    /* calculate radix-2 value equal or above window_length */
    window_order = 0;
    window_r2 = 1;
    while( window_r2 < fft_length) {
      window_order++;
      window_r2 <<= 1;
    }

    plhs[0] = MXCREATEFULL( window_r2, num_slices, COMPLEX);
    if (plhs[0] == NULL) {
      tfsaErr( "quadtfd", "Memory allocation failed: plhs[0]");
      winNTcheck( nlhs, plhs );
      return;
    }
    result_r = mxGetPr( plhs[0]);
    result_i = mxGetPi( plhs[0]);

    G_real = (double **)mxCalloc( size, sizeof(double *));
    if (G_real==NULL) {
      tfsaErr( "quadtfd", "Memory allocation failed: G_real");
      winNTcheck( nlhs, plhs );
      return;
    }
    G_real[0] = (double *)mxCalloc( (int)(size*size), sizeof(double));
    if (G_real[0]==NULL) {
      tfsaErr( "quadtfd", "Memory allocation failed");
      winNTcheck( nlhs, plhs );
      return;
    }
    for (i = 1; i < size; i++)
      G_real[i] = G_real[0] + (unsigned long)size * i;

    G_imag = (double **)mxCalloc( size, sizeof(double *));
    if (G_imag==NULL) {
      tfsaErr( "quadtfd", "Memory allocation failed: G_imag");
      winNTcheck( nlhs, plhs );
      return;
    }
    G_imag[0] = (double *)mxCalloc( (int)(size*size), sizeof(double));
    if (G_imag[0]==NULL) {
      tfsaErr( "quadtfd", "Memory allocation failed: G_imag[0]");
      winNTcheck( nlhs, plhs );
      return;
    }
    for (i = 1; i < size; i++)
      G_imag[i] = G_imag[0] + (unsigned long)size * i;

    complex_rm_kernel( G_real, window_length);

    /* G_imag is set to zero by Calloc above */


    quad_complexG_asymmG( signal_r, signal_i, signal_length, result_r, result_i,
			  num_slices, time_res,	window_length,
			  window_r2, window_order, G_real, G_imag);

    mxFree( (void *)(G_real[0]));
    mxFree( (void *)G_real);
    mxFree( (void *)(G_imag[0]));
    mxFree( (void *)G_imag);
  }
  else {      /* THIS IS DISABLED -- JR */
    /* distribution is defined by supplied kernel */
    int complex;
    int full_kernel;
    double *tmp_kernel;

    size = (unsigned)mxGetN( prhs[3]);
    full_kernel = (size == window_length);
    tmp_kernel = mxGetPr( prhs[3]);

    complex = mxIsComplex( prhs[3]);

    G_real = (double **)mxCalloc( size, sizeof(double *));
    if (G_real==NULL) {
      tfsaErr( "quadtfd", "Memory allocation failed: G_real");
      winNTcheck( nlhs, plhs );
      return;
    }
    G_real[0] = (double *)mxCalloc( (int)(size*size), sizeof(double));
    if (G_real[0]==NULL) {
      tfsaErr( "quadtfd", "Memory allocation failed: G_real[0]");
      winNTcheck( nlhs, plhs );
      return;
    }
    for (i = 1; i < size; i++)
      G_real[i] = G_real[0] + (unsigned long)size * i;

    for (i = 0; i < size; i++)
      for (j = 0; j< size; j++)
	G_real[i][j] = tmp_kernel[i + j * (unsigned long)size];

    if (complex) {
      tmp_kernel = mxGetPi( prhs[3]);
      G_imag = (double **)mxCalloc( size, sizeof(double *));
      if (G_imag==NULL) {
	tfsaErr( "quadtfd", "Memory allocation failed: G_imag");
	winNTcheck( nlhs, plhs );
	return;
      }
      G_imag[0] = (double *)mxCalloc( (int)(size*size), sizeof(double));
      if (G_imag[0]==NULL) {
	tfsaErr( "quadtfd", "Memory allocation failed: G_imag[0]");
	winNTcheck( nlhs, plhs );
	return;
      }
      for (i = 1; i < size; i++)
	G_imag[i] = G_imag[0] + (unsigned long)size * i;

      for (i = 0; i < size; i++)
	for (j = 0; j< size; j++)
	  G_imag[i][j] = tmp_kernel[i + j * (unsigned long)size];

      plhs[0] =	MXCREATEFULL( window_r2, num_slices, COMPLEX);
      if (plhs[0] == NULL) {
	tfsaErr( "quadtfd", "Memory allocation failed: plhs[0]");
	winNTcheck( nlhs, plhs );
	return;
      }
      result_r = mxGetPr( plhs[0]);
      result_i = mxGetPi( plhs[0]);

      if (full_kernel)
	quad_complexG_asymmG( signal_r, signal_i, signal_length,
			      result_r, result_i, num_slices, time_res,
			      window_length, window_r2, window_order,
			      G_real, G_imag);
      else
	quad_complexG_symmG( signal_r, signal_i, signal_length,
			     result_r, result_i, num_slices, time_res,
			     window_length, window_r2, window_order,
			     G_real, G_imag);

      mxFree( (void *)(G_imag[0]));
      mxFree( (void *)G_imag);
    }
    else {
      plhs[0] = MXCREATEFULL( window_r2, num_slices, REAL);
      if (plhs[0] == NULL) {
	tfsaErr( "quadtfd", "Memory allocation failed: plhs[0]");
	winNTcheck( nlhs, plhs );
	return;
      }
      result_r = mxGetPr( plhs[0]);

      if (full_kernel)
	quad_realG_asymmG( signal_r, signal_i, signal_length, result_r,
			   num_slices, time_res, window_length,
			   window_r2, window_order, G_real);
      else
	quad_realG_symmG( signal_r, signal_i, signal_length, result_r,
			  num_slices, time_res, window_length,
			  window_r2, window_order, G_real);
    }

  }
}




