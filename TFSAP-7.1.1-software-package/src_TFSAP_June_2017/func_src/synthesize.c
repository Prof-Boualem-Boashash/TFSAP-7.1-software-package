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
* Routine for synthesizing the signal from a given TFD
* Valid for WVD, STFT and Spectrogram.
*************************************************/

#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


#include "mex.h"


#include "arithm.h"		     /*	routines for complex data */
#include "tlocal.h"		     /*	local function prototypes */

#include "synthesize.h"
#include "analyt.h"
#include "fft.h"
#include "tfsa_c.h"
#include "window.h"
#include "spec.h"
#include "analyt.h"


#define	TFSA_PI	3.14159265358979323846

/* provide a fix for problem with DLL (windows) version of TFSA in the calculation of trigonometric functions.  For detailed explanation, see comments in sfpe.c */

#define	EPS 1e-5			/* arbitrary, but small	offset */

#ifdef DLLMEX				/* Add offset for PC version only... */
#define	COS(x)	cos( (double)(x	+ EPS) )
#define	SIN(x)	sin( (double)(x	+ EPS) )
#else
#define	COS(x) cos(x)
#define	SIN(x) sin(x)
#endif

/* Iteration limit for MSPEC algorithm ..*/
#define IT_LIMIT 500

/* Beginning of function */

int synthesize( int     analysis_type,
		double  *synth_signal_re,
		double  *synth_signal_im,
		int     signal_length,
		int     window_length,
		int     window_type,
		double  tol,
		double  *signal_real,
		double  *signal_imag,
		double  *tfd_re,
		double  *tfd_im,
		int     tfd_M,
		int     tfd_N,
		double  *rankone_err)
{
  unsigned i,l;
  int hlf_win;
  int smooth_win_width;
  double window_sum;
  double *window;
  double *y_real = NULL, *y_imag = NULL;
  double *X_re = NULL, *X_im = NULL;
  double *tfd_sqrt_re = NULL;
  double *A_re = NULL, *A_im = NULL;
  double *s_re = NULL, *s_im = NULL;
  double *u_re = NULL, *u_im = NULL;
  complex *synth_sig_even = NULL;
  double x_abs, err, norm_err, diag_mat11;
  int window_order, window_r2, M, M_spec, iterations;
  int a_M, hlf_sig, x, y;

  MATRIX *fig_hCancel, *fig_cancel, *fig_rhs[3], *plot_rhs[3], *fig_handle;
  char str_buffer[100], str_cancel[3];
  double hFigure, hCancel;

  /* TMP, db. */
/*   char *db_msg[50]; */


  /* Calculate radix 2 value */
  M = (int)pow( 2, ceil( log(tfd_M) / log(2) ) );


  /* Allocate memory. */
  window  = (double*) CALLOC( window_length, sizeof(double) );
  if( !window ){
    tfsaErr("synthesize", "Memory allocation failed: window");
    return 1;
  }
  y_real = (double*) CALLOC( M * tfd_N, sizeof(double) );
  if( !y_real ){
    tfsaErr("synthesize", "Memory allocation failed.");
    return 1;
  }
  y_imag = (double*) CALLOC( M * tfd_N, sizeof(double) );
  if( !y_imag ){
    tfsaErr("synthesize", "Memory allocation failed.");
    return 1;
  }




  /* Cal. Y = IFFT( tfd ) */
  if( analysis_type != MSPEC ){
    if( inverseFFT( tfd_re, tfd_im, tfd_M, tfd_N, M, y_real, y_imag ) )
      return 1;
  }



  /* Select and construct a smoothing window of the type specified by the user */
  smooth_win_width=window_length;
  hlf_win = ( window_length - 1 ) / 2;
  switch( window_type ) {

  case RECT:
    for (l=0;l<(unsigned)smooth_win_width;l++){
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
    for (l=1;l<(unsigned)smooth_win_width;l++) {
	window[l] = 1.0/(double)smooth_win_width;
    }
    break;
  }


  memset( synth_signal_re, 0, (size_t)( sizeof(double) * signal_length ) );
  memset( synth_signal_im, 0, (size_t)( sizeof(double) * signal_length ) );


  switch( analysis_type ){

    /* ------ */
    /*  STFT  */
    /* ------ */

  case IDFT:

    for( i=0; i<(unsigned)signal_length; i++ ){
      synth_signal_re[i] = *(y_real + i + i * M) / window[hlf_win];
      synth_signal_im[i] = *(y_imag + i + i * M) / window[hlf_win];
    }
    break;

  case OLA:

    window_sum = 0;
    for( i=0; i<(unsigned)window_length; i++ )
      window_sum += window[i];


    for( i=0; i<(unsigned)signal_length; i++){
      for( l=0; l<(unsigned)signal_length; l++){
	synth_signal_re[i] += *(y_real + i + l * M ) / window_sum;
	synth_signal_im[i] += *(y_imag + i + l * M ) / window_sum;
      }
    }
    break;

  case MSTFT:

    mstft( y_real, y_imag, M, synth_signal_re, synth_signal_im,
	   signal_length, window, window_length, hlf_win );

    break;


    /* ------------- */
    /*  SPECTROGRAM  */
    /* ------------- */

  case MSPEC:

          M_spec = M / 2 + 1;
          X_re  = (double *)CALLOC( M_spec * tfd_N, sizeof(double) );
          if( !X_re ) {
              tfsaErr("synthesize", "Memory allocation failed");
              return 1;
          }
          X_im  = (double *)CALLOC( M_spec * tfd_N, sizeof(double) );
          if( !X_im ) {
              tfsaErr("synthesize", "Memory allocation failed");
              return 1;
          }
          tfd_sqrt_re = (double *)CALLOC( tfd_M * tfd_N, sizeof(double) );
          if( !tfd_sqrt_re ) {
              tfsaErr("synthesize", "Memory allocation failed");
              return 1;
          }
          
          /* Initial estimate of signal */
          generateRandomSignal( synth_signal_re, synth_signal_im, signal_length );
          
          
          /* Calculate |S(t,f)| */
          norm_err = 0;
          for( i = 0; i < (unsigned)tfd_M; i++ ){
              for( l = 0; l < (unsigned)tfd_N; l++ ){
                  *(tfd_sqrt_re + i + l * tfd_M) = sqrt( (double)*(tfd_re + i + l * tfd_M) );
                  norm_err += (double)*(tfd_sqrt_re + i + l * tfd_M);
              }
          }
          
          /* FFT length (and order) needed for calculating STFT */
          window_order = 0;
          window_r2 = 1;
          while( window_r2 < M ){
              window_order ++;
              window_r2 <<= 1;
          }
          
          /* Set the error to some initial value */
          err = 2 * tol;
          
          
          /* Include an updating graph of error v. iterations */

          
          
          
          iterations = 0;
          
          while( (err > tol) && (iterations < IT_LIMIT) ) {
              
              iterations +=1;
              
              
              spec( synth_signal_re, synth_signal_im, signal_length, X_re, X_im,
                   signal_length, 1, window_length, window_type, window_r2,
                   window_order, (int)1 );
              
              
              err = 0;
              for( i = 0; i < (unsigned)M_spec; i++ ){
                  for( l = 0; l < (unsigned)tfd_N; l++ ){
                      
                      x_abs = module( cmplx(  *(X_re + i + l * M_spec) , *(X_im + i + l * M_spec) ) );
                      
                      if( x_abs != 0 ){
                          
                          *(X_re + i + l * M_spec) =  (double)*(tfd_sqrt_re + i + l * tfd_M) *
                          ( (double)*(X_re + i + l * M_spec) / x_abs );
                          *(X_im + i + l * M_spec) =  (double)*(tfd_sqrt_re + i + l * tfd_M) *
                          ( (double)*(X_im + i + l * M_spec) / x_abs );
                          
                      } else {
                          
                          *(X_re + i + l * M_spec) =  *(tfd_sqrt_re + i + l * tfd_M);
                          *(X_im + i + l * M_spec) =  0;
                          
                      }
                      
                      /* error is always one iteration behind but should suffice, */
                      /* saves computationally                                    */
                      err += ( x_abs - (double)*(tfd_sqrt_re + i + l * tfd_M) ) *
                      ( x_abs - (double)*(tfd_sqrt_re + i + l * tfd_M) );
                      
                  }
              }
              err = err / norm_err;
              
              
              inverseFFT( X_re, X_im, M_spec, tfd_N, M, y_real, y_imag );
              
              /* Need to clear the array here */
              memset( synth_signal_re, 0, (size_t)( sizeof(double) * signal_length ) );
              memset( synth_signal_im, 0, (size_t)( sizeof(double) * signal_length ) );
              
              
              mstft( y_real, y_imag, M, synth_signal_re, synth_signal_im, signal_length,
                    window, window_length, hlf_win );
              
              
              

              
          } /* end while( err > tol ) */
          

          
          /* If original signal is supplied then can correct phase on */
          /* synthesized signal.                                      */
          if( signal_real ){
              if( reconstructPhase( signal_real, signal_imag, signal_length,
                                   synth_signal_re, synth_signal_im ) )
                  return 1;
          }
          
          
          
          FREE((void *)X_re );
          FREE((void *)X_im );
          FREE((void *)tfd_sqrt_re);
    break;


    /* --------------------------- */
    /*  Wigner-Ville Distribution  */
    /* --------------------------- */

  case MWVD:

    /* Allocate the memory */
    a_M = (int)( tfd_N / 2 );
    hlf_sig = (int)( signal_length / 2 );
    A_re = (double*) CALLOC( a_M * a_M, sizeof(double) );
    if( !A_re ){
      tfsaErr("synthesize", "Memory allocation failed.");
      return 1;
    }
    A_im = (double*) CALLOC( a_M * a_M, sizeof(double) );
    if( !A_im ){
      tfsaErr("synthesize", "Memory allocation failed.");
      return 1;
    }
    s_re = (double*) CALLOC( a_M * a_M, sizeof(double) );
    if( !s_re ){
      tfsaErr("synthesize", "Memory allocation failed.");
      return 1;
    }
    s_im = (double*) CALLOC( a_M * a_M, sizeof(double) );
    if( !s_im ){
      tfsaErr("synthesize", "Memory allocation failed.");
      return 1;
    }
    u_re = (double*) CALLOC( a_M * a_M, sizeof(double) );
    if( !u_re ){
      tfsaErr("synthesize", "Memory allocation failed.");
      return 1;
    }
    u_im = (double*) CALLOC( a_M * a_M, sizeof(double) );
    if( !u_im ){
      tfsaErr("synthesize", "Memory allocation failed.");
      return 1;
    }
    synth_sig_even = (complex*) CALLOC( hlf_sig, sizeof(complex) );
    if( !synth_sig_even ){
      tfsaErr("synthesize", "Memory allocation failed.");
      return 1;
    }



    for( i = 0; i < (unsigned)a_M; i++ ){
      for( l = 0; l < (unsigned)a_M; l++ ){

	x = i - l;
	y = i + l;

	if( abs(x) <= window_length ){

	  /* Co-ords. in y are swapped as tfd matrix is represented that way */
	  if( x < 0 ){
	    x += tfd_M;
	    *(A_re + i + l * a_M) = *(y_real + x + y * tfd_M);
	    *(A_im + i + l * a_M) = *(y_imag + x + y * tfd_M);
	  } else {
	    *(A_re + i + l * a_M) = *(y_real + x + y * tfd_M);
	    *(A_im + i + l * a_M) = *(y_imag + x + y * tfd_M);
	  }
	}

      }
    }


    /* Singular value decomposition */
    if( do_svd( A_re, A_im, a_M, u_re, u_im, s_re, s_im ) )
      return 1;


    diag_mat11 = (double)*(s_re);
    if( ( diag_mat11 <= 0 ) || ( (s_im) && ((double)*(s_im)) ) ){
      tfsaErr("synthesize", "Matrix is not diagonalizable.");
      return 1;
    }
    diag_mat11 = sqrt( diag_mat11 );



    /* Synthesize the even samples from the signal */
    for( i = 0; i < (unsigned)hlf_sig; i++ ){
      synth_sig_even[i].re = diag_mat11 * (double)*(u_re + i);
      if( u_im )
	synth_sig_even[i].im = diag_mat11 * (double)*(u_im + i);
    }


    /* Calculate rank 1 approx. error */
    *rankone_err = 0;
    for( i = 1; i < (unsigned)a_M; i++ )
      *rankone_err += (double)*(s_re + i + i * a_M);



    /* Interpolate between the even samples */
    if( interpolate( synth_signal_re, synth_signal_im, signal_length, synth_sig_even ) )
      return 1;

    /* If original signal is supplied then can correct phase on */
    /* synthesized signal.                                      */
    if( signal_real ){
      if( reconstructPhase( signal_real, signal_imag, signal_length,
			    synth_signal_re, synth_signal_im ) )
	return 1;
    }



    /* Deallocate memory */
    FREE( (void *)A_re );
    FREE( (void *)A_im );
    FREE( (void *)u_re );
    FREE( (void *)s_re );
    if( u_im ) FREE( (void *)u_im );
    if( s_im ) FREE( (void *)s_im );
    if( synth_sig_even ) FREE( (void *)synth_sig_even );

    break;

  default:
    tfsa_cerr("Incorrect choice of signal synthesis routine");
    break;
  }



  /* Clear out memory */
  FREE((void *)window);
  FREE((void *)y_real);
  FREE((void *)y_imag);

  return 0;
}



/* --------------------------------------------------------------- */
/* Function that synthesizes the signal from a given modified STFT */
/* --------------------------------------------------------------- */

void
mstft(double *y_real,
      double *y_imag,
      int    M,
      double *synth_signal_re,
      double *synth_signal_im,
      int    signal_length,
      double *window,
      int    window_length,
      int    hlf_win)
{
  double window_energy;
  int i,l;


  window_energy = 0;
  for( i=0; i<window_length; i++ )
    window_energy += window[i] * window[i];


  for(i=0; i<signal_length; i++){
    for(l=0; l<signal_length; l++){

      if( (l - i + hlf_win >= 0) && ( l - i + hlf_win  < window_length ) ){
	synth_signal_re[i] += *(y_real + i + l * M ) * window[l - i + hlf_win];
	synth_signal_im[i] += *(y_imag + i + l * M ) * window[l - i + hlf_win];
      }
    }
    synth_signal_re[i] /= window_energy;
    synth_signal_im[i] /= window_energy;
  }

}



/* ---------------------------------------------------------- */
/* Calculates the IFFT of matrix (tfd_re, tfd_im) and outputs */
/* to (y_real, y_imag)                                        */
/* ---------------------------------------------------------- */

int
inverseFFT(double  *tfd_re,
	   double  *tfd_im,
	   int     tfd_M,
	   int     tfd_N,
	   int     M,
	   double  *y_real,
	   double  *y_imag)
{



  MATRIX *fft_lhs, *fft_rhs[2];



  if( tfd_im ){
    fft_rhs[0] = MXCREATEFULL (tfd_M, tfd_N, COMPLEX);
  } else {
    fft_rhs[0] = MXCREATEFULL (tfd_M, tfd_N, REAL);
  }
  fft_rhs[1] = MXCREATEFULL (1, 1, REAL);

  if ( !fft_rhs[0] || !fft_rhs[1] ){
    tfsaErr("synthesize", "Memory allocation failed" );
    return 1;
  }


  memcpy( (double *) mxGetPr( fft_rhs[0] ), tfd_re,
	  (size_t)( sizeof(double) * tfd_M * tfd_N ) );
  if( tfd_im ){
    memcpy( (double *) mxGetPi( fft_rhs[0] ), tfd_im,
	    (size_t)( sizeof(double) * tfd_M * tfd_N ) );
  }
  *( mxGetPr( fft_rhs[1] ) ) = M;


  /* Call FFT routine in MATLAB */
  mexCallMATLAB(1, &fft_lhs, 2, fft_rhs, "ifft" );
  mxDestroyArray( (MATRIX *)fft_rhs[0] );
  mxDestroyArray( (MATRIX *)fft_rhs[1] );


  memcpy( y_real , (double *)mxGetPr( fft_lhs ),
	  (size_t)( sizeof(double) * M * tfd_N ) );
  memcpy( y_imag , (double *)mxGetPi( fft_lhs ),
	  (size_t)( sizeof(double) * M * tfd_N ) );


  mxDestroyArray( (MATRIX *) fft_lhs );




  return 0;
}


/* ----------------------------------------- */
/* Generates a random signal between -1 -> 1 */
/* (only 6 point precision)                  */
/* ----------------------------------------- */

void
generateRandomSignal(double *sig_re,
		     double *sig_im,
		     int    signal_length)
{
  int i,j;
  time_t tm;
  tm = time( NULL );


  /* Set seed to changing variable (any int will do) */
  srand((int)tm);


  /* random number between -1 -> 1 with 6 point precision */
  for(i=0;i<signal_length;i++){
    sig_re[i] = 0;
    sig_im[i] = 0;

    for(j=1;j<6;j++){
      sig_re[i] += (double)((10.0*rand()/(RAND_MAX+1.0)) / pow(10,j) );
      sig_im[i] += (double)((10.0*rand()/(RAND_MAX+1.0)) / pow(10,j) );
    }

    if( (10.0*rand()/(RAND_MAX+1.0)) > 5 )
      sig_re[i] *= -1;

    if( (10.0*rand()/(RAND_MAX+1.0)) > 5 )
      sig_im[i] *= -1;

  }


}


/* ------------------------------------ */
/* Singular Value Decomposition routine */
/* Only returns U, S for USV = svd( A ) */
/* ------------------------------------ */

int
do_svd( double *A_re,
	double *A_im,
	int    a_M,
 	double *u_re,
	double *u_im,
	double *s_re,
	double *s_im)
{


  MATRIX *svd_lhs[3], *svd_rhs;


  svd_rhs = MXCREATEFULL (a_M, a_M, COMPLEX);

  if ( !svd_rhs ){
    tfsaErr("synthesize", "Memory allocation failed" );
    return 1;
  }


  memcpy( (double *) mxGetPr( svd_rhs ), A_re,
	  (size_t)( sizeof(double) * a_M * a_M ) );
  memcpy( (double *) mxGetPi( svd_rhs ), A_im,
	  (size_t)( sizeof(double) * a_M * a_M ) );



  /* Call Singular Value Decomposition routine in MATLAB */
  /* Returns matrices of the form : USV = svd( X )       */
  mexCallMATLAB(3, svd_lhs, 1, &svd_rhs, "svd" );
  mxDestroyArray( (MATRIX *)svd_rhs );


  memcpy( u_re , (double *)mxGetPr( svd_lhs[0] ),
	  (size_t)( sizeof(double) * a_M * a_M ) );

  if( mxGetPi( svd_lhs[0] ) ) {
    memcpy( u_im , (double *)mxGetPi( svd_lhs[0] ),
	    (size_t)( sizeof(double) * a_M * a_M ) );
  } else {
    u_im = NULL;
  }

  memcpy( s_re , (double *)mxGetPr( svd_lhs[1] ),
	  (size_t)( sizeof(double) * a_M * a_M ) );

  if( mxGetPi( svd_lhs[1] ) ){
    memcpy( s_im , (double *)mxGetPi( svd_lhs[1] ),
	    (size_t)( sizeof(double) * a_M * a_M ) );
  } else {
    s_im = NULL;
  }


  mxDestroyArray( (MATRIX *)svd_lhs[0] );
  mxDestroyArray( (MATRIX *)svd_lhs[1] );
  mxDestroyArray( (MATRIX *)svd_lhs[2] );


  return 0;
}


/* --------------------------------------------------------- */
/* Interpolation routine of order 2 for given complex signal */
/* --------------------------------------------------------- */

int
interpolate(double *synth_signal_re,
	    double *synth_signal_im,
	    int    signal_length,
	    complex *synth_sig_even)
{
  int i, iEven, S, S_order;
  complex *synth_sig_freq = NULL;


  synth_sig_freq = (complex*) CALLOC( signal_length, sizeof(complex) );
  if( !synth_sig_freq ){
    tfsaErr("synthesize", "Memory allocation failed.");
    return 1;
  }



  /* Reconstruct signal by interpolating */
  iEven = 0;
  for( i = 0; i < signal_length ; i++ ){
    if( !fmod(i, 2) ){
      synth_sig_freq[i] = synth_sig_even[iEven];
      iEven ++;
    }
  }

  /* signal_length should be radix 2 value already */
  S_order = (int) floor( log(signal_length) / log(2) );


  FFT( synth_sig_freq, S_order, 1 );

  /* In case signal_length was not radix 2 */
  S = (int)pow( 2, S_order );

  for( i = 0; i < S; i++ ){
    if( i < (S / 2) ){
      synth_sig_freq[i].re = synth_sig_freq[i].re * 2;
      synth_sig_freq[i].im = synth_sig_freq[i].im * 2;
    } else {
      synth_sig_freq[i] = cmplx( (double)0, (double)0 );
    }
  }

  FFT( synth_sig_freq, S_order, -1);

  for( i = 0; i < S; i++ ){
    synth_signal_re[i] = synth_sig_freq[i].re;
    synth_signal_im[i] = synth_sig_freq[i].im;
  }


  FREE( (void *)synth_sig_freq );

  return 0;
}


/* ------------------------------------------------------------ */
/* Calculates the phase from given signal and then reconstructs */
/* synthesized signal with this phase.                          */
/* ------------------------------------------------------------ */

int
reconstructPhase( double *original_sig_re,
		  double *original_sig_im,
		  int    signal_length,
		  double *synth_sig_re,
		  double *synth_sig_im)
{
  unsigned i;
  complex *sig = NULL, tmp, nom, denom, phase;


  /* Better to supply the analytic associate to avoid this.. */
  if( !original_sig_im ){
    sig = (complex *) CALLOC( signal_length, sizeof(complex) );
    if( !sig ) {
      tfsaErr("synthesize", "Memory allocation failed");
      return 1;
    }
    original_sig_im = (double *) CALLOC( signal_length, sizeof(double) );
    if( !original_sig_im ) {
      tfsaErr("synthesize", "Memory allocation failed");
      return 1;
    }



    for( i = 0; i < (unsigned)signal_length; i++ ){
      sig[i].re = original_sig_re[i];
    }

    default_sigana( sig, (unsigned)signal_length );


    for( i = 0; i < (unsigned)signal_length; i++ ){
      original_sig_im[i] = sig[i].im;
    }

  }



  nom = cmplx( (double)0, (double)0 );
  denom = cmplx( (double)0, (double)0 );
  tmp = cmplx( (double)0, (double)0 );


  /* Calculate the phase */
  for( i = 0; i < (unsigned)signal_length; i++ ){
    tmp = multpl( cmplx( original_sig_re[i], original_sig_im[i] ),
		  conj(cmplx( synth_sig_re[i] , synth_sig_im[i] )) );

    nom = add(  nom, tmp );
    denom = add( denom, cmplx( module( tmp ), 0 ) );
  }

  phase = divide( nom, denom );


  /* reconstruct signal with new phase */
  for( i = 0; i < (unsigned)signal_length; i++ ){
    tmp = multpl( cmplx( synth_sig_re[i], synth_sig_im[i] ) , phase );
    synth_sig_re[i] = tmp.re;
    synth_sig_im[i] = tmp.im;
  }


  if( sig ){
    FREE( (void *)sig );
    FREE( (void *)original_sig_im );
  }

  return 0;
}


