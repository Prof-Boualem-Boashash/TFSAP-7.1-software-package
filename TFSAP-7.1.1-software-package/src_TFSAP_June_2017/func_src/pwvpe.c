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
* Lastly the following reference is useful for understanding the basics of Instantaneous
* Frequency estimation:
* [4] B. Boashash, "Estimating and interpreting the instantaneous frequency of
* a signal—part 2: algorithms and applications", Proc. IEEE 80 (4) (1992) 540-568.
*
* Description:
*
* Routine for estimating the Instantaneous Frequency of a signal
* from	the peak of the	Polynomial Wigner-Ville	Distribution.
*************************************************/

#include <math.h>
#include <ctype.h>
#include <stdlib.h>

#include "mex.h"

#include "arithm.h"		     /*	routines for complex data */
#include "tlocal.h"		     /*	local function prototypes */

#include "pwvpe.h"
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

#define	TFSA_PI	3.14159265358979323846


#define	EPS 1e-5			/* arbitrary, but small	offset */

#ifdef DLLMEX							/* Add offset for PC version only... */
#define	COS(x)	cos( (double)(x	+ EPS) )
#define	SIN(x)	sin( (double)(x	+ EPS) )
#else
#define	COS(x) cos(x)
#define	SIN(x) sin(x)
#endif



void interpolate(int, int, complex*, int, int*, int*); 
void Maximum (complex*, complex*, int, int, double*, double*);

void pwvpe( double  *signal_real,	
           double  *signal_imag, 
           int     signal_length,
	   double  *result, 
           int     nplts, 
           int     time_res,
	   int     window_length, 
           int     window_r2,
	   int     window_order, 
           int	   deg)
{
  int i, j, j1,	k, max,	signal_length_r2, sliding_time;
  complex *sig,	*Wig, Kernel, X, Y, temp;
  double sampling_freq,	posmax;
  int hlf;
  double maximum;
  double lag, frac;
  int order, deg_order,	i1, i2,	i3, i4;
  int new_window_r2, new_window_order, new_window_length;
  int new_signal_length, new_time_res, new_length;


  /* Calculate radix 2 above or	equal signal length */
  order	= 0;
  signal_length_r2 = 1;
  while	( signal_length_r2 < signal_length) {
     signal_length_r2 <<= 1;
     order++;
  }

  /* Allocate complex array storage for	the input signal */
  sig =	(complex *) mxCalloc (signal_length_r2*deg, sizeof(complex));
  if (sig==NULL)
  {
    tfsaErr( "pwvpe", "Memory allocation failed");
    return;
  }

  Wig =	(complex *) mxCalloc (window_r2, sizeof(complex));
  if (Wig == NULL)
  {
    tfsaErr( "pwvpe", "Memory allocation failed");
    return;
  }


  /* Read in input signal into complex array */
  if (signal_imag == NULL) {
     for (i = 0; i < signal_length; i++) {
        sig[i].re = signal_real[i];
        sig[i].im = 0.0;
     }
    /* Now generate analytic signal of the real	signal */
    default_sigana(sig,	signal_length);
  } else {
    for	(i=0; i	< signal_length; i++) {
      sig[i].re	= signal_real[i];
      sig[i].im	= signal_imag[i];
    }
  }

  /*
   * Now perform interpolation on the signal in	the time domain
   * to	get the	nearest	to the time lag	values of the 4th order
   * Kernel (0.794)
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
  new_length = (int)pow(2.,(double)order);
   
 
  /*
   * Begin calculating PWVD for	two windows simulataneously
   */
  frac = (double)(1.0/pow(2.0,(double)(1.0/3.0)));
  sliding_time = 0;
  for (i=0; i< nplts -1; i+=2) {  

    hlf	= (int) new_window_length/2;
    if ( (sliding_time - hlf) < 0) {
      hlf = sliding_time;
    }
    if (sliding_time+hlf >= new_signal_length) {
      hlf = new_signal_length -1 -sliding_time;
    }

    for	(j= ((int) hlf/deg) + 1; j < window_r2 - (int) hlf/deg; j++) {
      Wig[j].re	= 0.0;
      Wig[j].im	= 0.0;
    }

    /* form kernel by computing	mutiplied shifted versions of
     * the signal at all values	of time	lag for	every time slide
     */
    for	(j=0; j	<= hlf ; j+=deg) {

      /* Indices for two Kernel	terms that occur at integer lags */
      i1 = sliding_time	+ j;
      i2 = sliding_time	- j;

      /* Indices for two Kernel	terms that occur at fraction lags */
      lag = (double)(frac * (double)j);
      i3 = sliding_time	+ (int)lag;
      i4 = sliding_time	- (int)lag;

      /* Mutiplication of integer lag terms */
      X.re = sig[i1].re	* sig[i2].re + sig[i1].im * sig[i2].im;
      X.im =-sig[i2].re	* sig[i1].im + sig[i1].re * sig[i2].im;

      /* Multiplication	and squaring of	noninteger lag terms */
      Y.re = sig[i3].re	* sig[i4].re + sig[i3].im * sig[i4].im;
      Y.im = sig[i3].im	* sig[i4].re - sig[i3].re * sig[i4].im;
      temp.re =	Y.re * Y.re - Y.im * Y.im;
      temp.im =	Y.re * Y.im * 2.0;

      /* Now form the Kernel */
      Kernel.re	= X.re * temp.re - X.im	* temp.im;
      Kernel.im	= X.re * temp.im + X.im	* temp.re;

      /* Prepare Kernel	for FFT	*/
      j1 = (int) j/deg;
      Wig[j1].re = Kernel.re;
      Wig[j1].im = Kernel.im;
      if (j != 0) {
	Wig[window_r2-j1].re = Kernel.re;
	Wig[window_r2-j1].im = -Kernel.im;
      }
    }

    if (i < nplts-1) {
      sliding_time += new_time_res;

      /* now adjust effective time window length to ensure that	window is
       *  contained within signal
       */
      hlf = (int) new_window_length/2;
      if ( sliding_time - hlf < 0) {
	hlf = sliding_time;
      }

      if ((sliding_time + hlf) >= new_signal_length) {
	hlf = new_signal_length	-1 -sliding_time;
      }

      /* form kernel by	computing mutiplied shifted versions of
       * the signal at all values of time lag for every	time slide
       */
      for (j=0;	j <= hlf ; j+=deg)
      {
	/* Indices for two Kernel terms	that occur at integer lags */
	i1 = sliding_time + j;
	i2 = sliding_time - j;

	/*Indices for two Kernel terms that occur at fraction lags */
	lag = (double)(frac * (double)j);
	i3 = sliding_time + (int)lag;
	i4 = sliding_time - (int)lag;

	/* Mutiplication of integer lag	terms */
	X.re = sig[i1].re * sig[i2].re + sig[i1].im * sig[i2].im;
	X.im =-sig[i2].re * sig[i1].im + sig[i1].re * sig[i2].im;

	/* Multiplication and squaring of noninteger lag terms */
	Y.re = sig[i3].re * sig[i4].re + sig[i3].im * sig[i4].im;
	Y.im = sig[i3].im * sig[i4].re - sig[i3].re * sig[i4].im;
	temp.re	= Y.re * Y.re -	Y.im * Y.im;
	temp.im	= Y.re * Y.im *	(double)2;

	/* Now form the	Kernel */
	Kernel.re = X.re * temp.re - X.im * temp.im;
	Kernel.im = X.re * temp.im + X.im * temp.re;

	/* Prepare Kernel for FFT */
	j1 =  (int) j/deg;
	Wig[j1].re += -Kernel.im;
	Wig[j1].im +=  Kernel.re;
	if (j1 != 0)
	{
	  Wig[window_r2-j1].re += Kernel.im;
	  Wig[window_r2-j1].im += Kernel.re;
	}
      }
    }

    FFT	(Wig, window_order, 1);

    /*
     * Find out	the frequency value at which PWIG
     * is Maximum. Before we have to frequency scale
     * the resulting distribution
     */
    max	= (int)(window_r2*(double)(1.0/0.85)*0.5);

    sampling_freq = 1.;

    /* Find out	the
     * frequency value at which	PWIG is	Maximum.
     */
    maximum = Wig[0].re;

    posmax =0.;
    for	(k=1; k < window_r2; k++)
      if (Wig[k].re > maximum) {
	  maximum = Wig[k].re;
	  posmax = (double) k;
      }
    *(result + i) = sampling_freq * posmax / ((1.0/0.85)*(double)window_r2);

    /* Find out	the
     * frequency value at which	PWIG is	Maximum.
     */
    maximum = Wig[0].im;
    posmax =0.;
    for	(k=1; k<window_r2; k++) {
      if (Wig[k].im > maximum) {
	  maximum = Wig[k].im;
	  posmax = (double)k;
      }
    }
    *(result + i+1) = sampling_freq * posmax / ((1.0/0.85)*(double)window_r2);
    
    sliding_time += new_time_res;
  }

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
    sig[i].re =	(double)0;
    sig[i].im =	(double)0;
  }

  /* Fast Fourier Transform the	zero padded analytic signal */
  FFT (sig, *order, 1);

  for ( i = (int) signal_length_r2/2;  i < signal_length_r2; i++)
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

