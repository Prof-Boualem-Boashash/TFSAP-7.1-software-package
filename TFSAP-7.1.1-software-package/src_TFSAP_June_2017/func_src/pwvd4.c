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
* signal characteristics, SoftwareX, 2017s.
* In addition, the following 3rd reference is relevant for use of the executable
* [3] B. Boashash (ed.), Time-Frequency Signal Analysis and Processing, 2nd Edition
* (London: Elsevier / Academic Press, December 2015); ISBN 978-0-12-398499-9. Book website:
* http://booksite.elsevier.com/9780123984999/
*
* Description:
*
* Routine for generating fourth order Kernel Polynomial Wigner-Ville
* Distribution
*
* This	program	is a highly  optimised version,	taking account of the
* shape of the	kernel for the pwvd, and the reality of	the resulting
* distribution.
*************************************************/

#include <math.h>
#include <ctype.h>
#include <stdlib.h>

#include "mex.h"

#include "arithm.h"		     /*	routines for complex data */
#include "tlocal.h"		     /*	local function prototypes */

#include "pwvd4.h"
#include "analyt.h"
#include "fft.h"
#include "tfsa_c.h"

#ifdef WAITBAR
MATRIX *prhs[2];
MATRIX *p1, *p2;
MATRIX *plhs[1];
double *n1;
double waitbar_hndl;
#endif


int pwvd4( double *signal_real,	double *signal_imag, int
	signal_length, double *result, int nplts, int
	time_res, int window_length,
	int window_r2, int window_order)
{
  int i, j, max, deg, order, signal_length_r2;
  complex *sig,	*Wig, X, Y, Kernel;
  int new_window_r2, new_window_length,	new_signal_length;
  int new_time_res;
  int new_window_order;
  int sliding_time, hlf, i1, i2;
  void interpolate(int,	int, complex*, int,int*);


  /* Degree of interpolation.
   */
  deg =	2;

  /* Calculate radix 2 above or	equal signal length */
  order	= 0;
  signal_length_r2 = 1;
  while	( signal_length_r2 < signal_length)
  {
	signal_length_r2 <<= 1;
	order++;
  }

  /* Allocate complex array storage for	the input signal */
  max =	signal_length_r2*deg;
  sig =	(complex *) mxCalloc (max, sizeof(complex));
  if (sig==NULL)
    {
      tfsaErr("pwvd", "Memory allocation failed" );
      return 1;
    }
  /* Read in input signal into complex array */
  if (signal_imag == NULL)
  {
    for	(i = 0;	i < signal_length; i++)
    {
      sig[i].re	= signal_real[i];
      sig[i].im	= 0.0;
    }

    /* Now generate analytic signal of the real	signal */
    default_sigana(sig,signal_length);
  }
  else
    for	(i=0; i	< signal_length; i++)
    {
      sig[i].re	= signal_real[i];
      sig[i].im	= signal_imag[i];
    }

  /*
   * Now perform interpolation on the signal in	the time domain
   * to	prevent	the ocurrence of aliasing problem.
   */
  interpolate(signal_length_r2,	signal_length, sig, deg, &order);

  /*
   * Calculate new parameters for new scaling
   */
  new_window_r2	= deg *	window_r2;
  new_window_order = window_order+1; /*	1=log2(deg) */
  new_window_length = window_length;
  new_signal_length = deg * signal_length;
  new_time_res = deg * time_res;


  /* Allocate storage for both windows of PWVD */
  max =	new_window_r2;
  Wig	= (complex *) mxCalloc (max, sizeof(complex));
  if (Wig==NULL)
  {
     tfsaErr( "pwvd", "Memory allocation failed" );
     return 1;
  }

  /*
   * Begin calculating PWVD for	two windows simulataneously
   */

#ifdef WAITBAR
  p1 = mxCreateFull(1, 1, REAL);
  p2 = mxCreateString("Please wait, computing distribution...");
  n1 = mxGetPr(p1);
  *n1 =	0.0;
  prhs[0] = p1;
  prhs[1] = p2;
  mexCallMATLAB( 1, plhs, 2, prhs, "waitbar");
  waitbar_hndl = *mxGetPr(plhs[0]);
#endif

  sliding_time = 0;
  for (i=0; i<nplts; i+=2)
  {  /*
      *	Computes the Polynomial	Wigner-Ville distribution of two windows of
      *	the analytic signal sig	simultaneously.	 This is done because it
      *	is known that a	WVD is real and	hence the complex part of the
      *	array Wig passed to FFT	can be used to compute the second window.
      *
      *
      *
      *	Note: Complex Array Wig	contains the kernel of window 1	in its
      *	      real part, and the kernel	of window 2 in its imaginary
      *	      part, in order to	perform	one FFT	rather than two. This
      *	      saving is	possible because both kernel and its DFT are real.
      *
      */

#ifdef WAITBAR
  *n1 =	(double)i/nplts;
  mexCallMATLAB( 0, plhs, 1, prhs, "waitbar");
#endif

    /* NOTE: we	do two frames at a time.  The second frame is
     *	multiplied by j	and added to the first.	 The resulting spectrum
     *	after FFT contains the FFT of the first	frame and second frame
     *	in its real and	imaginary parts	respectively.
     *
     * now adjust effective time window	length to ensure that window is
     *	contained within signal
     */
    hlf	= new_window_length/2;
    if ((int)sliding_time - (int)hlf < 0)
    {
      hlf = sliding_time;
    }
    if (sliding_time+hlf >= new_signal_length)
    {
      hlf = new_signal_length -1 -sliding_time;
    }

    for	(j=hlf+1; j < new_window_r2 - hlf; j++)
    {
      Wig[j].re	= 0.0;
      Wig[j].im	= 0.0;
    }

    /* form kernel by computing	mutiplied shifted versions of
     * the signal at all values	of time	lag for	every time slide
     */
    for	(j=0; j	<= hlf ; j++)
    {
      /* Indices for two Kernel	terms that occur at integer lags */
      i1 = sliding_time	+ j;
      i2 = sliding_time	- j;

      /* Squaring and mutiplication of integer lag terms */
      X.re = sig[i1].re	* sig[i1].re - sig[i1].im * sig[i1].im;
      X.im = sig[i1].re	* sig[i1].im * (float)2;
      Y.re = sig[i2].re	* sig[i2].re - sig[i2].im * sig[i2].im;
      Y.im =-sig[i2].im	* sig[i2].re * (float)2;

      /* Now form the Kernel */
      Kernel.re	= X.re * Y.re -	X.im * Y.im;
      Kernel.im	= X.re * Y.im +	X.im * Y.re;

      /* Prepare Kernel	for FFT	*/
      Wig[j].re	= Kernel.re;
      Wig[j].im	= Kernel.im;
      if (j != 0)
      {
	Wig[new_window_r2-j].re	= Kernel.re;
	Wig[new_window_r2-j].im	=-Kernel.im;
      }
    }

    if (i < nplts-1)
    {
      sliding_time += new_time_res;

      /* now adjust effective time window length to ensure that	window is
       *  contained within signal
       */
      hlf = new_window_length/2;
      if ((int)sliding_time - (int)hlf < 0)
      {
	hlf = sliding_time;
      }
      if (sliding_time+hlf >= new_signal_length)
      {
	hlf = new_signal_length	-1 -sliding_time;
      }

      /* form kernel by	computing mutiplied shifted versions of
       * the signal at all values of time lag for every	time slide
       */
      for (j=0;	j <= hlf ; j++)
      {
	/* Indices for two Kernel terms	that occur at integer lags */
	i1 = sliding_time + j;
	i2 = sliding_time - j;

	/* Squaring and	mutiplication of integer lag terms */
	X.re = sig[i1].re * sig[i1].re - sig[i1].im * sig[i1].im;
	X.im = sig[i1].re * sig[i1].im * (float)2;
	Y.re = sig[i2].re * sig[i2].re - sig[i2].im * sig[i2].im;
	Y.im =-sig[i2].im * sig[i2].re * (float)2;

	/* Now form the	Kernel */
	Kernel.re = X.re * Y.re	- X.im * Y.im;
	Kernel.im = X.re * Y.im	+ X.im * Y.re;

	/* Prepare Kernel for FFT */
	Wig[j].re +=-Kernel.im;
	Wig[j].im += Kernel.re;
	if (j != 0)
	{
	  Wig[new_window_r2-j].re += Kernel.im;
	  Wig[new_window_r2-j].im += Kernel.re;
	}
      }
    }

    FFT(Wig, new_window_order, 1);

    /* Now store the first slice of WVT	after correcting the magnitude
     * of the spectrum with a factor of	2, because of the interpolation
     * see Oppenheim & Schafer,	Discrete-Time Signal Processing,
     * sec 3.6.2
     */
    for(j=0; j < (int)(window_r2*(float)4*0.5);	j++)
    {
      *(result + i*(int)(window_r2*(float)4*0.5)+j) = 16.*Wig[j].re;
				      /* 16=pow(2,4) is	interpolation scale */
      if (i < nplts-1)
	 *(result+(i+1)*(int)(window_r2*(float)4*0.5)+j) = 16.*Wig[j].im;
    }

    sliding_time += new_time_res;
  }

#ifdef WAITBAR
  *n1 =	(double)1;
  mexCallMATLAB( 0, plhs, 1, prhs, "waitbar");
  *n1 =	waitbar_hndl;
  mexCallMATLAB( 0, plhs, 1, prhs, "close");
  mxFreeMatrix(p1);
  mxFreeMatrix(p2);
#endif


  return 0;
}



void interpolate(int signal_length_r2, int signal_length,
complex	*sig, int deg, int *order)

/*
 * Interpolating the signal in the time	domain
 */

  {
  int alpha = 0, beta =	1;
  int max, i;

  /* degree of interpolation must be radix 2 */
  while	(beta <	deg)
  {
    beta <<= 1;
    alpha++;
  }
  deg =	beta;
  max =deg * signal_length_r2;
  /* Zero padding */
  for (i=signal_length;	i<max; i++)
  {
    sig[i].re =	(double)0;
    sig[i].im =	(double)0;
  }

  /* Fast Fourier Transform the	zero padded analytic signal */
  FFT (sig, *order, 1);

  for (i=signal_length_r2/2; i<signal_length_r2; i++)
  {
    sig[i+(deg-1)*signal_length_r2].re = sig[i].re;
    sig[i+(deg-1)*signal_length_r2].im = sig[i].im;
    sig[i].re =	(double)0;
    sig[i].im =	(double)0;
  }
  *order += alpha;
  FFT (sig, *order, -1);
}







