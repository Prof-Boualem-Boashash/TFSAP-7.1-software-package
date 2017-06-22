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
* Routine for calculating the spectrogram representation of a signal.
* Implemented as the instantaneous windowed periodogram estimate
*
* Notes: This code was taken from the function "sfpe"
* (peak of the spectrogram IF estimate) and adapted to provide a
* spectrogram using the same form of input variables as the
* "quadratic" function.
*************************************************/

#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>


#include "mex.h"


#include "arithm.h"		     /*	routines for complex data */
#include "tlocal.h"		     /*	local function prototypes */

#include "spec.h"
#include "analyt.h"
#include "fft.h"
#include "tfsa_c.h"
#include "window.h"

#define	TFSA_PI	3.14159265358979323846

/* provide a fix for problem with DLL (windows) version of TFSA in the calculation of trigonometric functions. */

#define	EPS 1e-5			/* arbitrary, but small	offset */

#ifdef DLLMEX				/* Add offset for PC version only... */
#define	COS(x)	cos( (double)(x	+ EPS) )
#define	SIN(x)	sin( (double)(x	+ EPS) )
#else
#define	COS(x) cos(x)
#define	SIN(x) sin(x)
#endif

/* Beginning of function */

int spec( double  *signal_real, 
	    double  *signal_imag, 
	    int	    signal_length,
	    double  *tfd_re,
	    double  *tfd_im,
	    int     numslices, 
	    int     time_res,
	    int     window_length, 
	    int     window_type,
	    int     window_r2,
	    int     window_order,
	    unsigned stft_or_spec)
{
  unsigned i, k, l, n;
  int m, M;
  int hlf;
  int smooth_win_width;
  int z_index;
  int sl;
  complex *sig;
  double *window;
  double *time_lag_kernel_re = NULL, *time_lag_kernel_im = NULL;


  MATRIX *fft_lhs, *fft_rhs;
  double *stft_re = NULL, *stft_im = NULL;



  /* Allocate complex array storage for	the input signal */
  sig =	(complex *) CALLOC(window_r2, sizeof(complex));
  if( !sig ) {
    tfsaErr( "spec", "Memory allocation failed");
    return 1;
  }


  /* Allocate memory for the window */
  window  = (double*) CALLOC(window_length, sizeof(double));
  if( !window ) {
    tfsaErr("spec", "Memory allocation failed: window");
    return 1;
  }


  /* Allocate memory for the signal kernel */
  time_lag_kernel_re  = (double *)CALLOC( window_r2 * numslices, sizeof(double) );
  if( !time_lag_kernel_re ) {
    tfsaErr("spec", "Memory allocation failed: window");
    return 1;
  }
  time_lag_kernel_im  = (double *)CALLOC( window_r2 * numslices, sizeof(double) );
  if( !time_lag_kernel_im ) {
    tfsaErr("spec", "Memory allocation failed: window");
    return 1;
  }
  stft_re = (double *)CALLOC( window_r2 * numslices, sizeof(double) );
  if( !stft_re ) {
    tfsaErr("spec", "Memory allocation failed");
    return 1;
  }
  stft_im = (double *)CALLOC( window_r2 * numslices, sizeof(double) );
  if( !stft_im ) {
    tfsaErr("spec", "Memory allocation failed");
    return 1;
  }



  /* Signal is extended periodically in window_r2
   * when cal. STFT for each window. 
   * Therefore need to truncate or pad signal from signal_length
   */
  for (i = 0; i < window_r2; i++) {
      sig[i].re = 0;
      sig[i].im = 0;
  }

  /* Find the min value between length and period */
  sl = ( window_r2 < signal_length ? window_r2 : signal_length );

  /* Read in input signal into complex array */
  if (signal_imag == NULL) {
    
    /* initialise the complex holding array */ 
    for (i = 0; i < sl; i++) {
	  sig[i].re = signal_real[i];
    }
      
    /* Now generate analytic signal of the real signal by */
    /* sending the signal out to analyt function.         */
    default_sigana(sig, window_r2);  
      
  } else {
  
    for (i=0; i < sl; i++) {
      sig[i].re = signal_real[i];
      sig[i].im = signal_imag[i];
    }

  }


  /* Select and construct a smoothing window of the type specified by the user */
  smooth_win_width=window_length;
  switch( window_type ) {

  case RECT:
    for (l=0;l<(int)smooth_win_width;l++){
      window[l] = 1.0/(double)smooth_win_width;
    }
    break;
  case HANN:
    hann(window, smooth_win_width);
    break;
  case HAMM:
    hamming(window, smooth_win_width);
    break;
  case BART:
    bartlett(window, smooth_win_width);
    break;
  default:
    for (l=1;l<(int)smooth_win_width;l++) {
      window[l] = 1.0/(double)smooth_win_width;
    }
    break;

  }

  
  /* Calculate Spectogram */


  hlf = window_length/2; 
  n = 0;
  for(k = 0; k < numslices; k++){
    
    /* Form time-lag kernel */
    for (m = -hlf; m <= hlf; m++){
      
 	  /* z[.] periodic in 'window_r2' */
      z_index = fmod( (window_r2 + m + n), window_r2 );
      
      
      *(time_lag_kernel_re + (z_index) + (k * window_r2)) = 
	sig[z_index].re * window[m+hlf];
      *(time_lag_kernel_im + (z_index) + (k * window_r2)) = 
	sig[z_index].im * window[m+hlf];
      
    }
    n += time_res; 
  }


  /* Use MATLAB's FFT function as quicker than C code for matrices */
  fft_rhs = MXCREATEFULL( window_r2, numslices, COMPLEX );
  if ( !fft_rhs ){ 
    tfsaErr("spec", "Memory allocation failed" );
    return 1;
  }


  memcpy( (double *)mxGetPr( fft_rhs ), time_lag_kernel_re, 
	  (size_t)( sizeof(double) * window_r2 * numslices ) );
  memcpy( (double *)mxGetPi( fft_rhs ), time_lag_kernel_im, 
	  (size_t)( sizeof(double) * window_r2 * numslices ) );


  /* Call FFT routine in MATLAB */
  if( mexCallMATLAB(1, &fft_lhs, 1, &fft_rhs, "fft" ) ){
    tfsaErr("spec", "Unable to call FFT MATLAB function" );
    return 1;
  }

  mxDestroyArray( (MATRIX *)fft_rhs ); 


  memcpy( stft_re , (double *)mxGetPr( fft_lhs ),
	  (size_t)( sizeof(double) * window_r2 * numslices ) );
  memcpy( stft_im , (double *)mxGetPi( fft_lhs ),
	  (size_t)( sizeof(double) * window_r2 * numslices ) );

  mxDestroyArray( (MATRIX *)fft_lhs ); 



  M = (int)( window_r2 / 2 + 1 );

  /* If calcuating the Spectrogram then take the |STFT|^2 */
  if( !stft_or_spec ){

    for( n = 0; n < M; n++ ){
      for( m = 0; m < numslices; m++ ){
	*(tfd_re + n + m * M)  = ( *(stft_re + n + m * window_r2) )  * 
	  ( *(stft_re + n + m * window_r2) ) +
	  ( *(stft_im + n + m * window_r2) ) *
	  ( *(stft_im + n + m * window_r2) );
      }
    }

  } else {

    for( n = 0; n < M; n++ ){
      for( m = 0; m < numslices; m++ ){
	*(tfd_re + n + m * M)  = ( *(stft_re + n + m * window_r2) );
	*(tfd_im + n + m * M)  = ( *(stft_im + n + m * window_r2) );
      }
    }

  }



  /* Free memory */
  FREE((void *)sig);
  FREE((void *)window);
  FREE((void *)stft_re);
  FREE((void *)stft_im);
  /* The MATLAB function will free this memory */
  /* as the memory was passed to a MATLAB matrix */
  FREE((void *)time_lag_kernel_re);
  FREE((void *)time_lag_kernel_im);







  return 0;
}

