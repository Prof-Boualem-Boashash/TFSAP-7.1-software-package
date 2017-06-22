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
* Routine for calculating the Rihaczek and Levin distributions directly.
*************************************************/

#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>


#include "mex.h"


#include "arithm.h"		     /*	routines for complex data */
#include "tlocal.h"		     /*	local function prototypes */

#include "rihaczek.h"
#include "analyt.h"
#include "fft.h"
#include "tfsa_c.h"
#include "window.h"
#include "spec.h"

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

int rihaczek( double  *signal_real, 
	      double  *signal_imag, 
	      int     signal_length,
	      double  *result_re,
	      double  *result_im,
	      int     numslices, 
	      int     time_res,
	      int     signal_r2,
	      int     signal_order,
	      int     type_of_dist,
	      int     window_length,
	      int     window_type)
{
  unsigned i, k, n, m;
  unsigned sl, hlf_sig_len;
  complex *sig, *ft_sig, result_tmp;
  double exp_tmp_re, exp_tmp_im, exp_tmp_arg;
  double *stft_re = NULL, *stft_im = NULL;


  /* Allocate complex array storage for	the input signal */

  sig =	(complex *) CALLOC( signal_r2, sizeof(complex) );
  if ( !sig ){ 
    tfsaErr( "rihaczek", "Memory allocation failed");
    return 1;
  }


  /* Allocate complex array storage for	the input signal */
  ft_sig = (complex *) CALLOC( signal_r2, sizeof(complex) );
  if( !ft_sig ){
    tfsaErr( "rihaczek", "Memory allocation failed");
    return 1;
  }



  /* Signal is extended periodically in signal_r2 .
   * Therefore need to truncate or pad signal from signal_length
   */
  for (i = 0; i < signal_r2; i++) {
      sig[i].re = 0;
      sig[i].im = 0;
  }

  /* Find the min value between length and period */
  sl = ( signal_r2 < signal_length ? signal_r2 : signal_length );

  /* Read in input signal into complex array */
  if (signal_imag == NULL) 
    {
      for (i = 0; i < sl; i++) 
	{
	  /* initialise the complex holding array */ 
	  sig[i].re = signal_real[i];
	}
      
      /* Now generate analytic signal of the real signal by */
      /* sending the signal out to analyt function.         */
      default_sigana(sig, signal_r2);  
      
    } 
  else 
    { 
      for (i=0; i < sl; i++) 
	{
	  sig[i].re = signal_real[i];
	  sig[i].im = signal_imag[i];
	}
    }

  hlf_sig_len = ( signal_r2 / 2) + 1; 

  /* Clear out TFD matrix */
  for(k = 0; k < numslices; k++) {
    for(m = 0; m < hlf_sig_len; m++) {
      *(result_re + m + k*hlf_sig_len) = 0;
      if( result_im )
	*(result_im + m + k*hlf_sig_len) = 0;
    }
  }


  


  /* Calculate windowed - Rihaczek Distribution */
  if( window_length > 0 ){

    /* First need to calculate (STFT(z(t)))^* */
    stft_re  = (double *)CALLOC( hlf_sig_len * numslices, sizeof(double) );
    if( !stft_re ) {
      tfsaErr("rihaczek", "Memory allocation failed");
      return 1;
    }
    stft_im  = (double *)CALLOC( hlf_sig_len * numslices, sizeof(double) );
    if( !stft_im ) {
      tfsaErr("rihaczek", "Memory allocation failed");
      return 1;
    }

    spec( signal_real, signal_imag, signal_length, stft_re, stft_im, 
	  numslices, time_res, window_length, window_type, signal_r2,
	  signal_order, (int)1 );


    n = 0;
    for(i = 0; i < numslices; i++) {

      
      for (k = 0; k < hlf_sig_len; k++){

	exp_tmp_arg = ( (2 * TFSA_PI * n * k) / signal_r2 );
	exp_tmp_re = cos( exp_tmp_arg ); 
	exp_tmp_im = -sin( exp_tmp_arg ); 
	
	/* R[n,k] = z[n]S^*[n,k].exp(-j.2.pi.n.m/N), where S[n,k] if STFT of z[n] */
	
	result_tmp =  
	  multpl( (multpl( sig[n], conj( cmplx( *(stft_re + k + (i * hlf_sig_len)) , *(stft_im + k + (i * hlf_sig_len)) ) ) )),  (cmplx( exp_tmp_re, exp_tmp_im )) );  
	
	/* Fill in TFD matrix */
	*(result_re + k + (i * hlf_sig_len)) =  result_tmp.re;
	if( result_im )
	  *(result_im + k + (i * hlf_sig_len)) =  result_tmp.im;
      }
      
      n += time_res; 
    }

    FREE( (void *)stft_re );
    FREE( (void *)stft_im );

  /* Calculate Rihaczek Distribution */
  } else {

 
    for(i=0; i< signal_r2; i++){
      ft_sig[i] = sig[i];
    }

    /* Find Z(f) of z(t) */
    FFT( ft_sig, signal_order, 1 );

    n = 0;
    for(i = 0; i < numslices; i++) 
      {
      
	for (k = 0; k < hlf_sig_len; k++)
	  {
	    exp_tmp_arg = ( (2 * TFSA_PI * n * k) / signal_r2 );
	    exp_tmp_re = cos( exp_tmp_arg ); 
	    exp_tmp_im = -sin( exp_tmp_arg ); 
	    result_tmp =  
	      multpl( (multpl( sig[n], conj( ft_sig[k] ) )), (cmplx( exp_tmp_re, exp_tmp_im )) );  

	    /* Fill in TFD matrix */
	    *(result_re + k + (i * hlf_sig_len)) =  result_tmp.re;
	    if( result_im )
	      *(result_im + k + (i * hlf_sig_len)) =  result_tmp.im;
	  }

	n += time_res; 
      }
  }


  /* Clear out memory */
  FREE((void *)sig);
  FREE((void *)ft_sig);
  
  return 0;
}

