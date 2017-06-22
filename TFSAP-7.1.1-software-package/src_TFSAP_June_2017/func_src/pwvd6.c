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
* Routine for generating Sixth	order Kernel Polynomial	Wigner-Ville
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

#include "pwvd6.h"
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

int pwvd61( double *signal_real,	double *signal_imag, int
	signal_length, double *result, int nplts, int
	time_res, int window_length, int window_r2,
	int window_order, int deg)
{
  int i, j, j1,	max, order, signal_length_r2;
  complex *sig,	*Wig, X1, X2, X3, Kernel;
  int new_window_r2, new_window_length,	new_signal_length;
  int new_time_res, new_length;
  int new_window_order,	deg_order;
  int sliding_time, hlf, i1, i2, i3, i4, i5, i6;
  float	lag1, lag2, lag3, frac1, frac2, frac3;
  void interpolate(int,	int, complex*, int, int*, int*);


  /* Degree of interpolation.
   * note: any change of this value must take into account
   *	   the size of the output array	(result) in the
   *	   Gateway routine pwvd_m.c
   *
   */

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
    tfsaErr( "pwvd61", "Memory allocation failed");
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
   * Initialize	Waitbar
   */

#ifdef WAITBAR
  p1 = mxCreateFull(1, 1, REAL);
  p2 = mxCreateString("Please wait, Computing Distribution...");
  n1 = mxGetPr(p1);
  *n1 =	0.0;
  prhs[0] = p1;
  prhs[1] = p2;
  mexCallMATLAB( 1, plhs, 2, prhs, "waitbar");
  waitbar_hndl = *mxGetPr(plhs[0]);
#endif

  /*
   * Now perform interpolation on the signal in	the time domain
   * to	get the	nearest	to the time lag	values of the 6th order
   * Kernel.
   */
  interpolate(signal_length_r2,	signal_length, sig, deg, &order, &deg_order);

  /*
   * Calculate new parameters for new scaling
   */
  new_window_r2	= deg *	window_r2;
  new_window_order = window_order+deg_order;
  new_window_length = deg * window_length;
  new_signal_length = deg * signal_length;
  new_time_res = deg * time_res;
  new_length = (int)pow(2.,(float)order);

  /* Allocate storage for both windows of PWVD */
  max =	window_r2;
  Wig	= (complex *) mxCalloc (max, sizeof(complex));
  if (Wig==NULL)
  {
   tfsaErr( "pwvd61", "Memory allocation	failed");
   return 1;
  }

  /*
   * Begin calculating PWVD for	two windows simulataneously
   */
    frac1 = 0.62;
    frac2 = 0.75;
    frac3 = 0.87;
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

     /*
      *	NOTE: we do two	frames at a time.  The second frame is
      *	 multiplied by j and added to the first.  The resulting	spectrum
      *	 after FFT contains the	FFT of the first frame and second frame
      *	 in its	real and imaginary parts respectively.
      *
      *	now adjust effective time window length	to ensure that window is
      *	 contained within signal
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

    for	(j=(hlf/deg)+1;	j < window_r2 -	(hlf/deg); j++)
    {
      Wig[j].re	= 0.0;
      Wig[j].im	= 0.0;
    }

    /* form kernel by computing	mutiplied shifted versions of
     * the signal at all values	of time	lag for	every time slide
     */
    for	(j=0; j	<= hlf ; j+=deg)
    {
      /* Indices for two Kernel	terms that occur at fraction lags (0.62)*/
      lag1 = (float)(frac1 * (float)j);
      i1 = sliding_time	+ (int)lag1;
      i2 = sliding_time	- (int)lag1;

      /* Indices for two Kernel	terms that occur at fraction lags (0.75)*/
      lag2 = (float)(frac2 * (float)j);
      i3 = sliding_time	+ (int)lag2;
      i4 = sliding_time	- (int)lag2; 

      /* Indices for two Kernel	terms that occur at fraction lags (-0.87)*/
      lag3 = (float)(frac3 * (float)j);
      i5 = sliding_time	- (int)lag3;
      i6 = sliding_time	+ (int)lag3; 

      /* Mutiplication of lag1 terms */
      X1.re = sig[i1].re  * sig[i2].re + sig[i1].im * sig[i2].im;
      X1.im =-sig[i1].re  * sig[i2].im + sig[i1].im * sig[i2].re;

      /* Mutiplication of lag2 terms */
      X2.re = sig[i3].re  * sig[i4].re + sig[i3].im * sig[i4].im;
      X2.im =-sig[i3].re  * sig[i4].im + sig[i3].im * sig[i4].re;

      /* Mutiplication of lag3 terms */
      X3.re = sig[i5].re  * sig[i6].re + sig[i5].im * sig[i6].im;
      X3.im =-sig[i5].re  * sig[i6].im + sig[i5].im * sig[i6].re;

      /* Now form the Kernel */
      Kernel.re	= X1.re * X2.re * X3.re - X1.im * X2.im * X3.re - X1.re * X2.im * X3.im - X1.im * X2.re * X3.im;
      Kernel.im	= X1.re * X2.re * X3.im - X1.im * X2.im * X3.im + X1.re * X2.im * X3.re + X1.im * X2.re * X3.re;

      /* Prepare Kernel	for FFT	*/
      j1=j/deg;
      Wig[j1].re = Kernel.re;
      Wig[j1].im = Kernel.im;
      if (j != 0)
      {
	Wig[window_r2-j1].re = Kernel.re;
	Wig[window_r2-j1].im = -Kernel.im;
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
      for (j=0;	j <= hlf ; j+=deg)
      {
	/* Indices for two Kernel terms	that occur at fractional lags (0.62) */
	lag1 = (float)(frac1 * (float)j);
	i1 = sliding_time + (int)lag1;
	i2 = sliding_time - (int)lag1;

	/* Indices for two Kernel terms	that occur at fractional lags (0.75) */
	lag2 = (float)(frac2 * (float)j);
	i3 = sliding_time + (int)lag2;
	i4 = sliding_time - (int)lag2;

	/* Indices for two Kernel terms	that occur at fractional lags (-0.87) */
	lag3 = (float)(frac3 * (float)j);
	i5 = sliding_time - (int)lag3;
	i6 = sliding_time + (int)lag3;

      /* Mutiplication of lag1 terms */
      X1.re = sig[i1].re  * sig[i2].re + sig[i1].im * sig[i2].im;
      X1.im =-sig[i1].re  * sig[i2].im + sig[i1].im * sig[i2].re;

      /* Mutiplication of lag2 terms */
      X2.re = sig[i3].re  * sig[i4].re + sig[i3].im * sig[i4].im;
      X2.im =-sig[i3].re  * sig[i4].im + sig[i3].im * sig[i4].re;

      /* Mutiplication of lag3 terms */
      X3.re = sig[i5].re  * sig[i6].re + sig[i5].im * sig[i6].im;
      X3.im =-sig[i5].re  * sig[i6].im + sig[i5].im * sig[i6].re;


	/* Now form the	Kernel */
	Kernel.re = X1.re * X2.re * X3.re - X1.im * X2.im * X3.re - X1.re * X2.im * X3.im - X1.im * X2.re * X3.im;
	Kernel.im = X1.re * X2.re * X3.im - X1.im * X2.im * X3.im + X1.re * X2.im * X3.re + X1.im * X2.re * X3.re;

	/* Prepare Kernel for FFT */
	j1=j/deg;
	Wig[j1].re +=-Kernel.im;
	Wig[j1].im += Kernel.re;
	if (j1 != 0)
	{
	  Wig[window_r2-j1].re += Kernel.im;
	  Wig[window_r2-j1].im += Kernel.re;
	}
      }
    }

    FFT(Wig, window_order, 1);

    for(j=0; j < (int)(window_r2*0.5); j++)
    {
      *(result + i*(int)(window_r2*0.5)+j) = Wig[j].re ;
      if (i < nplts-1)
	 *(result+(i+1)*(int)(window_r2*0.5)+j) = Wig[j].im ;
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

  mxFree((void *)sig);
  mxFree((void *)Wig);

  return 0;
}


void interpolate(int signal_length_r2, int signal_length,complex *sig,
				 int deg, int *order,int *dgordr)
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
    sig[i].re =	(double)0.0;
    sig[i].im =	(double)0.0;
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
  *dgordr = alpha;
  FFT (sig, *order, -1);
}







