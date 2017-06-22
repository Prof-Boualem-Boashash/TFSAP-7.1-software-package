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
*
* quadtfd
*
* Quadratic TF Distributions
*
*************************************************/

#include <stdlib.h>
#include <math.h>

#include "mex.h"

#include "arithm.h"		     /*	routines for complex data */
#include "tlocal.h"		     /*	local function prototypes */

#include "quadtfd.h"
#include "fft.h"
#include "analyt.h"
#include "tfsa_c.h"

int Modulus_Index(int , int );
void analytic_signal( double *signal_real,
        double *z_real_out, double *z_imag_out, int signal_length);

/* Global variables for	waitbar	*/
#ifdef WAITBAR
static MATRIX *prhs[2];
static MATRIX *p1, *p2;
static MATRIX *plhs[1];
static double *n1;
static double waitbar_hndl;
#endif

int quad_complexG_asymmG( double     *signal_real,
                           double     *signal_imag,
	                   unsigned   signal_length,
                           double     *result_real,
                           double     *result_imag,
                           unsigned   num_slices,
                           unsigned   time_res,
	                   unsigned   window_length,
                           unsigned   window_r2,
                           unsigned   window_order,
                           double     **G_real,
                           double     **G_imag)
{
  unsigned i, j, n;
  complex *sig,	*wvd;
  unsigned w_centre;
  unsigned hlf;
  unsigned i1, i2, i3, i4;
  double X, Y;

  /* read in input signal into complex array */
  sig =	(complex *) mxCalloc (signal_length, sizeof(complex));
  if (sig==NULL) {
    tfsaErr( "quadtfd","Memory allocation	failed"	);
    return 1;
  }

  if (signal_imag == NULL) {
    for	(i = 0;	i < signal_length; i++)	{
      sig[i].re	= signal_real[i];
      sig[i].im	= 0.0;
    }
    default_sigana(sig,signal_length);
  }
  else
    for	(i=0; i	< signal_length; i++) {
      sig[i].re	= signal_real[i];
      sig[i].im	= signal_imag[i];
    }

  /* allocate storage for slice	of wvd */
  wvd =	(complex *) mxCalloc (window_r2, sizeof(complex));

  if (wvd==NULL) {
    tfsaErr( "quadtfd", "Memory allocation failed" );
    mxFree((void *)sig);
    return 1;
  }

  /* starting at 0 is not real clever, since the first slice will
     apply a window of 1 long.	For want of a better procedure,	and by
     the decision of the powers	that be, we will start at zero.	*/

  w_centre = 0;

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

  for (i = 0; i	< num_slices; i++) {

#ifdef WAITBAR
    *n1	= (double)i/num_slices;
    mexCallMATLAB( 0, plhs, 1, prhs, "waitbar");
#endif

    /* now adjust effective window length to ensure that window	is
       contained within	signal */

    /* note that window	length across slices (for 2d kernel) will be
       the same	as the window length along kernel */

    hlf	= window_length/2;
    if ((int)w_centre-(int)hlf < 0) {
      hlf = w_centre;
      /* effective_wl =	w_centre * 2 + 1; */
    }
    if (w_centre+hlf>=signal_length) {
      hlf = signal_length-1-w_centre;
      /* effective_wl =	(signal_length-1-w_centre) * 2 + 1; */
    }

    /* initialize wvd array (only those	points not explicitly set
       below */
    for	(j=hlf+1; j < window_r2-hlf; j++) {
      wvd[j].re	= 0.0;
      wvd[j].im	= 0.0;
    }

    /* do j = 0	loop separately	*/

    i1 = w_centre;
    i2 = w_centre;
    X =	sig[i1].re * sig[i2].re	+ sig[i1].im * sig[i2].im;
    Y =	sig[i2].re * sig[i1].im	- sig[i1].re * sig[i2].im;

    wvd[0].re =	G_real[0][0] * X - G_imag[0][0]	* Y;
    wvd[0].im =	G_real[0][0] * Y + G_imag[0][0]	* X;

    /* Now do:
       wvd[j] += G( n, j) ** ( sig[i1] ** conj(sig[i2])	+ sig[-i2] ** conj(sig[-i1]));
       */

    for	(n=1; n<=(hlf);	n++) {

      i1 = w_centre + n;
      i2 = w_centre + n;
      i3 = w_centre - n;
      i4 = w_centre - n;

      X	= sig[i1].re * sig[i2].re + sig[i1].im * sig[i2].im;
      Y	= sig[i2].re * sig[i1].im - sig[i1].re * sig[i2].im;

      wvd[0].re	+= G_real[window_length-n][0] *	X - G_imag[window_length-n][0] * Y;
      wvd[0].im	+= G_real[window_length-n][0] *	Y + G_imag[window_length-n][0] * X;

      X	= sig[i3].re * sig[i4].re + sig[i3].im * sig[i4].im;
      Y	= sig[i4].re * sig[i3].im - sig[i3].re * sig[i4].im;

      wvd[0].re	+= G_real[n][0]	* X - G_imag[n][0] * Y;
      wvd[0].im	+= G_real[n][0]	* Y + G_imag[n][0] * X;

    }

    /* do for all other	lags */

    /* form kernel */
    for	(j=1; j	<= hlf ; j++) {

      i1 = w_centre + j;
      i2 = w_centre - j;
      X	= sig[i1].re * sig[i2].re + sig[i1].im * sig[i2].im;
      Y	= sig[i2].re * sig[i1].im - sig[i1].re * sig[i2].im;

      wvd[j].re	= G_real[0][j] * X - G_imag[0][j] * Y;
      wvd[j].im	= G_real[0][j] * Y + G_imag[0][j] * X;

      wvd[window_r2-j].re = G_real[0][window_length-j] * X - G_imag[0][window_length-j]	* -Y;
      wvd[window_r2-j].im = G_real[0][window_length-j] * -Y + G_imag[0][window_length-j] * X;

      /* Now do:
      wvd[j] +=	G( n, j) ** ( sig[i1] ** conj(sig[i2]) + sig[-i2] ** conj(sig[-i1]));
      */

      for (n=1;	n<=(hlf-j); n++) {

	i1 = w_centre +	n + j;
	i2 = w_centre +	n - j;
	i3 = w_centre -	n + j;
	i4 = w_centre -	n - j;

	X = sig[i1].re * sig[i2].re + sig[i1].im * sig[i2].im;
	Y = sig[i2].re * sig[i1].im - sig[i1].re * sig[i2].im;

	wvd[j].re += G_real[window_length-n][j]	* X - G_imag[window_length-n][j] * Y;
	wvd[j].im += G_real[window_length-n][j]	* Y + G_imag[window_length-n][j] * X;

	wvd[window_r2-j].re += G_real[window_length-n][window_length-j]	* X -
	  G_imag[window_length-n][window_length-j] * -Y;
	wvd[window_r2-j].im += G_real[window_length-n][window_length-j]	* -Y +
	  G_imag[window_length-n][window_length-j] * X;

	X = sig[i3].re * sig[i4].re + sig[i3].im * sig[i4].im;
	Y = sig[i4].re * sig[i3].im - sig[i3].re * sig[i4].im;

	wvd[j].re += G_real[n][j] * X -	G_imag[n][j] * Y;
	wvd[j].im += G_real[n][j] * Y +	G_imag[n][j] * X;

	wvd[window_r2-j].re += G_real[n][window_length-j] * X -	G_imag[n][window_length-j] * -Y;
	wvd[window_r2-j].im += G_real[n][window_length-j] * -Y + G_imag[n][window_length-j] * X;
      }
    }

    FFT(wvd, window_order, 1);

    if (result_imag == NULL)
      for(j=0; j < window_r2; j++)
	*(result_real+i*window_r2+j) =	wvd[j].re;
    else
      for(j=0; j < window_r2; j++) {
	*(result_real+i*window_r2+j) = wvd[j].re;
	*(result_imag+i*window_r2+j) = wvd[j].im;
      }

    w_centre +=	time_res;
  }

#ifdef WAITBAR
  *n1 =	(double)1;
  mexCallMATLAB( 0, plhs, 1, prhs, "waitbar");
 *n1 = waitbar_hndl;
  mexCallMATLAB( 0, plhs, 1, prhs, "close");
  mxFreeMatrix(p1);
  mxFreeMatrix(p2);
#endif

  mxFree((void *)sig);
  mxFree((void *)wvd);

  return 0;
}


int quad_realG_asymmG(	double *signal_real,
                        double *signal_imag,
	                unsigned signal_length,
                        double	*result_real,
                        unsigned num_slices,
                        unsigned time_res,
                        unsigned window_length,
                        unsigned window_r2,
                        unsigned window_order,
                        double **G_real)
{
  unsigned i, j, n;
  complex *sig,	*wvd;
  unsigned w_centre;
  unsigned hlf;
  unsigned i1, i2, i3, i4;
  double X, Y;

  /* read in input signal into complex array */
  sig =	(complex *) mxCalloc (signal_length, sizeof(complex));
  if (sig==NULL) {
    tfsaErr( "quadtfd", "Memory allocation failed" );
    return 1;
  }

  if (signal_imag == NULL) {
    for	(i = 0;	i < signal_length; i++)	{
      sig[i].re	= signal_real[i];
      sig[i].im	= 0.0;
    }
    default_sigana(sig,signal_length);
  }
  else
    for	(i=0; i	< signal_length; i++) {
      sig[i].re	= signal_real[i];
      sig[i].im	= signal_imag[i];
    }

  /* allocate storage for slice	of wvd */
  wvd =	(complex *) mxCalloc (window_r2, sizeof(complex));
  if (wvd==NULL) {
    tfsaErr( "quadtfd","Memory allocation	failed"	);
    mxFree((void *)sig);
    return 1;
  }

  /* starting at 0 is not real clever, since the first slice will
     apply a window of 1 long.	For want of a better procedure,	and by
     the decision of the powers	that be, we will start at zero.	*/

  w_centre = 0;

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

  for (i = 0; i	< num_slices; i+=2) {

#ifdef WAITBAR
    *n1	= (double)i/num_slices;
    mexCallMATLAB( 0, plhs, 1, prhs, "waitbar");
#endif

    /* now adjust effective window length to ensure that window	is
       contained within	signal */

    /* note that window	length across slices (for 2d kernel) will be
       the same	as the window length along kernel */

    hlf	= window_length/2;
    if ((int)w_centre-(int)hlf < 0) {
      hlf = w_centre;
      /* effective_wl =	w_centre * 2 + 1; */
    }
    if (w_centre+hlf>=signal_length) {
      hlf = signal_length-1-w_centre;
      /* effective_wl =	(signal_length-1-w_centre) * 2 + 1; */
    }

    /* initialize wvd array (only those	points not explicitly set
       below */
    for	(j=hlf+1; j < window_r2-hlf; j++) {
      wvd[j].re	= 0.0;
      wvd[j].im	= 0.0;
    }

    /* do j = 0	loop separately	*/

    i1 = w_centre;
    i2 = w_centre;
    X =	sig[i1].re * sig[i2].re	+ sig[i1].im * sig[i2].im;
    Y =	sig[i2].re * sig[i1].im	- sig[i1].re * sig[i2].im;

    wvd[0].re =	G_real[0][0] * X;
    wvd[0].im =	G_real[0][0] * Y;

    /* Now do:
       wvd[j] += G( n, j) ** ( sig[i1] ** conj(sig[i2])	+ sig[-i2] ** conj(sig[-i1]));
       */

    for	(n=1; n<=(hlf);	n++) {

      i1 = w_centre + n;
      i2 = w_centre + n;
      i3 = w_centre - n;
      i4 = w_centre - n;

      X	= sig[i1].re * sig[i2].re + sig[i1].im * sig[i2].im;
      Y	= sig[i2].re * sig[i1].im - sig[i1].re * sig[i2].im;

      wvd[0].re	+= G_real[window_length-n][0] *	X;
      wvd[0].im	+= G_real[window_length-n][0] *	Y;

      X	= sig[i3].re * sig[i4].re + sig[i3].im * sig[i4].im;
      Y	= sig[i4].re * sig[i3].im - sig[i3].re * sig[i4].im;

      wvd[0].re	+= G_real[n][0]	* X;
      wvd[0].im	+= G_real[n][0]	* Y;

    }

    /* do for all other	lags */

    /* form kernel */
    for	(j=1; j	<= hlf ; j++) {

      i1 = w_centre + j;
      i2 = w_centre - j;
      X	= sig[i1].re * sig[i2].re + sig[i1].im * sig[i2].im;
      Y	= sig[i2].re * sig[i1].im - sig[i1].re * sig[i2].im;

      wvd[j].re	= G_real[0][j] * X;
      wvd[j].im	= G_real[0][j] * Y;

      wvd[window_r2-j].re = G_real[0][window_length-j] * X;
      wvd[window_r2-j].im = G_real[0][window_length-j] * -Y;

      /* Now do:
      wvd[j] +=	G( n, j) ** ( sig[i1] ** conj(sig[i2]) + sig[-i2] ** conj(sig[-i1]));
      */

      for (n=1;	n<=(hlf-j); n++) {

	i1 = w_centre +	n + j;
	i2 = w_centre +	n - j;
	i3 = w_centre -	n + j;
	i4 = w_centre -	n - j;

	X = sig[i1].re * sig[i2].re + sig[i1].im * sig[i2].im;
	Y = sig[i2].re * sig[i1].im - sig[i1].re * sig[i2].im;

	wvd[j].re += G_real[window_length-n][j]	* X;
	wvd[j].im += G_real[window_length-n][j]	* Y;

	wvd[window_r2-j].re += G_real[window_length-n][window_length-j]	* X;
	wvd[window_r2-j].im += G_real[window_length-n][window_length-j]	* -Y;

	X = sig[i3].re * sig[i4].re + sig[i3].im * sig[i4].im;
	Y = sig[i4].re * sig[i3].im - sig[i3].re * sig[i4].im;

	wvd[j].re += G_real[n][j] * X;
	wvd[j].im += G_real[n][j] * Y;

	wvd[window_r2-j].re += G_real[n][window_length-j] * X;
	wvd[window_r2-j].im += G_real[n][window_length-j] * -Y;
      }
    }

    if (i < num_slices-1) {

      w_centre += time_res;

      /* now adjust effective window length to ensure that window is
	 contained within signal */

      /* note that window length across	slices (for 2d kernel) will be
	 the same as the window	length along kernel */

      hlf = window_length/2;
      if ((int)w_centre-(int)hlf < 0) {
	hlf = w_centre;
	/* effective_wl	= w_centre * 2 + 1; */
      }
      if (w_centre+hlf>=signal_length) {
	hlf = signal_length-1-w_centre;
	/* effective_wl	= (signal_length-1-w_centre) * 2 + 1; */
      }

      /* do j =	0 loop separately */

      i1 = w_centre;
      i2 = w_centre;
      X	= sig[i1].re * sig[i2].re + sig[i1].im * sig[i2].im;
      Y	= sig[i2].re * sig[i1].im - sig[i1].re * sig[i2].im;

      wvd[0].re	+= G_real[0][0]	* -Y;
      wvd[0].im	+= G_real[0][0]	* X;

      /* Now do:
	 wvd[j]	+= G( n, j) ** ( sig[i1] ** conj(sig[i2]) + sig[-i2] **	conj(sig[-i1]));
	 */

      for (n=1;	n<=(hlf); n++) {

	i1 = w_centre +	n;
	i2 = w_centre +	n;
	i3 = w_centre -	n;
	i4 = w_centre -	n;

	X = sig[i1].re * sig[i2].re + sig[i1].im * sig[i2].im;
	Y = sig[i2].re * sig[i1].im - sig[i1].re * sig[i2].im;

	wvd[0].re += G_real[window_length-n][0]	* -Y;
	wvd[0].im += G_real[window_length-n][0]	* X;

	X = sig[i3].re * sig[i4].re + sig[i3].im * sig[i4].im;
	Y = sig[i4].re * sig[i3].im - sig[i3].re * sig[i4].im;

	wvd[0].re += G_real[n][0] * -Y;
	wvd[0].im += G_real[n][0] * X;

      }

      /* do for	all other lags */

      /* form kernel */
      for (j=1;	j <= hlf ; j++)	{

	i1 = w_centre +	j;
	i2 = w_centre -	j;
	X = sig[i1].re * sig[i2].re + sig[i1].im * sig[i2].im;
	Y = sig[i2].re * sig[i1].im - sig[i1].re * sig[i2].im;

	wvd[j].re += G_real[0][j] * -Y;
	wvd[j].im += G_real[0][j] * X;

	wvd[window_r2-j].re += G_real[0][window_length-j] * Y;
	wvd[window_r2-j].im += G_real[0][window_length-j] * X;

	/* Now do:
	   wvd[j] += G(	n, j) ** ( sig[i1] ** conj(sig[i2]) + sig[-i2] ** conj(sig[-i1]));
	   */

	for (n=1; n<=(hlf-j); n++) {

	  i1 = w_centre	+ n + j;
	  i2 = w_centre	+ n - j;
	  i3 = w_centre	- n + j;
	  i4 = w_centre	- n - j;

	  X = sig[i1].re * sig[i2].re +	sig[i1].im * sig[i2].im;
	  Y = sig[i2].re * sig[i1].im -	sig[i1].re * sig[i2].im;

	  wvd[j].re += G_real[window_length-n][j] * -Y;
	  wvd[j].im += G_real[window_length-n][j] * X;

	  wvd[window_r2-j].re += G_real[window_length-n][window_length-j] * Y;
	  wvd[window_r2-j].im += G_real[window_length-n][window_length-j] * X;

	  X = sig[i3].re * sig[i4].re +	sig[i3].im * sig[i4].im;
	  Y = sig[i4].re * sig[i3].im -	sig[i3].re * sig[i4].im;

	  wvd[j].re += G_real[n][j] * -Y;
	  wvd[j].im += G_real[n][j] * X;

	  wvd[window_r2-j].re += G_real[n][window_length-j] * Y;
	  wvd[window_r2-j].im += G_real[n][window_length-j] * X;
	}
      }
    }

    FFT(wvd, window_order, 1);

    for(j=0; j < window_r2; j++) {
      *(result_real+i*window_r2+j) =  wvd[j].re;
      if (i < num_slices-1)
	*(result_real+(i+1)*window_r2+j) = wvd[j].im;
    }

    w_centre +=	time_res;
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
  mxFree((void *)wvd);

  return 0;
}


int quad_complexG_symmG( double     *signal_real,
                          double     *signal_imag,
	                  unsigned   signal_length,
                          double     *result_real,
                          double     *result_imag,
                          unsigned   num_slices,
                          unsigned   time_res,
	                  unsigned   window_length,
                          unsigned   window_r2,
                          unsigned   window_order,
                          double     **G_real,
                          double     **G_imag)
{
  unsigned i, j, n;
  complex *sig,	*wvd;
  unsigned w_centre;
  unsigned hlf;
  unsigned i1, i2, i3, i4;
  double X, Y;

  /* read in input signal into complex array */
  sig =	(complex *) mxCalloc (signal_length, sizeof(complex));
  if (sig==NULL) {
    tfsaErr( "quadtfd", "Memory allocation failed");
    return 1;
  }

  if (signal_imag == NULL) {
    for	(i = 0;	i < signal_length; i++)	{
      sig[i].re	= signal_real[i];
      sig[i].im	= 0.0;
    }
    default_sigana(sig,signal_length);
  }
  else
    for	(i=0; i	< signal_length; i++) {
      sig[i].re	= signal_real[i];
      sig[i].im	= signal_imag[i];
    }

  /* allocate storage for slice	of wvd */
  wvd =	(complex *) mxCalloc (window_r2, sizeof(complex));
  if (wvd==NULL) {
    tfsaErr( "quadtfd","Memory allocation	failed");
    mxFree((void *)sig);
    return 1;
  }

  /* starting at 0 is not real clever, since the first slice will
     apply a window of 1 long.	For want of a better procedure,	and by
     the decision of the powers	that be, we will start at zero.	*/

  w_centre = 0;

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

  for (i = 0; i	< num_slices; i++) {

#ifdef WAITBAR
      *n1 = (double)i/num_slices;
      mexCallMATLAB( 0,	plhs, 1, prhs, "waitbar");
#endif

    /* now adjust effective window length to ensure that window	is
       contained within	signal */

    /* note that window	length across slices (for 2d kernel) will be
       the same	as the window length along kernel */

    hlf	= window_length/2;
    if ((int)w_centre-(int)hlf < 0) {
      hlf = w_centre;
      /* effective_wl =	w_centre * 2 + 1; */
    }
    if (w_centre+hlf>=signal_length) {
      hlf = signal_length-1-w_centre;
      /* effective_wl =	(signal_length-1-w_centre) * 2 + 1; */
    }

    /* initialize wvd array (only those	points not explicitly set
       below */
    for	(j=hlf+1; j < window_r2-hlf; j++) {
      wvd[j].re	= 0.0;
      wvd[j].im	= 0.0;
    }

    /* do j = 0	loop separately	*/

    i1 = w_centre;
    i2 = w_centre;
    X =	sig[i1].re * sig[i2].re	+ sig[i1].im * sig[i2].im;
    Y =	sig[i2].re * sig[i1].im	- sig[i1].re * sig[i2].im;

    wvd[0].re =	G_real[0][0] * X - G_imag[0][0]	* Y;
    wvd[0].im =	G_real[0][0] * Y + G_imag[0][0]	* X;

    /* Now do:
       wvd[j] += G( n, j) ** ( sig[i1] ** conj(sig[i2])	+ sig[-i2] ** conj(sig[-i1]));
       */

    for	(n=1; n<=(hlf);	n++) {

      i1 = w_centre + n;
      i2 = w_centre + n;
      i3 = w_centre - n;
      i4 = w_centre - n;

      X	= sig[i1].re * sig[i2].re + sig[i1].im * sig[i2].im;
      Y	= sig[i2].re * sig[i1].im - sig[i1].re * sig[i2].im;
      X	+= sig[i3].re *	sig[i4].re + sig[i3].im	* sig[i4].im;
      Y	+= sig[i4].re *	sig[i3].im - sig[i3].re	* sig[i4].im;

      wvd[0].re	+= G_real[n][0]	* X - G_imag[n][0] * Y;
      wvd[0].im	+= G_real[n][0]	* Y + G_imag[n][0] * X;

    }

    /* do for all other	lags */

    /* form kernel */
    for	(j=1; j	<= hlf ; j++) {

      i1 = w_centre + j;
      i2 = w_centre - j;
      X	= sig[i1].re * sig[i2].re + sig[i1].im * sig[i2].im;
      Y	= sig[i2].re * sig[i1].im - sig[i1].re * sig[i2].im;

      wvd[j].re	= G_real[0][j] * X - G_imag[0][j] * Y;
      wvd[j].im	= G_real[0][j] * Y + G_imag[0][j] * X;

      /* Now do:
      wvd[j] +=	G( n, j) ** ( sig[i1] ** conj(sig[i2]) + sig[-i2] ** conj(sig[-i1]));
      */

      for (n=1;	n<=(hlf-j); n++) {

	i1 = w_centre +	n + j;
	i2 = w_centre +	n - j;
	i3 = w_centre -	n + j;
	i4 = w_centre -	n - j;

	X = sig[i1].re * sig[i2].re + sig[i1].im * sig[i2].im;
	Y = sig[i2].re * sig[i1].im - sig[i1].re * sig[i2].im;
	X += sig[i3].re	* sig[i4].re + sig[i3].im * sig[i4].im;
	Y += sig[i4].re	* sig[i3].im - sig[i3].re * sig[i4].im;

	wvd[j].re += G_real[n][j] * X -	G_imag[n][j] * Y;
	wvd[j].im += G_real[n][j] * Y +	G_imag[n][j] * X;
      }

      wvd[window_r2-j].re = wvd[j].re;
      wvd[window_r2-j].im = -wvd[j].im;
    }

    FFT(wvd, window_order, 1);

    if (result_imag == NULL)
      for(j=0; j < window_r2; j++)
	*(result_real+i*window_r2+j) =	wvd[j].re;
    else
      for(j=0; j < window_r2; j++) {
	*(result_real+i*window_r2+j) = wvd[j].re;
	*(result_imag+i*window_r2+j) = wvd[j].im;
      }

    w_centre +=	time_res;
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
  mxFree((void *)wvd);

  return 0;
}


int quad_realG_symmG( double     *signal_real,
                       double     *signal_imag,
	               unsigned   signal_length,
                       double	  *result_real,
                       unsigned	  num_slices,
                       unsigned   time_res,
                       unsigned   window_length,
	               unsigned   window_r2,
                       unsigned   window_order,
                       double     **G_real)
{

    int N, i, j, f1, f2, t1, a, b, p, g1, g2, g3, g4, c, d, ind;
    double *z_real, *z_imag, *data1, *data2, *data3, *data4;

	int size_window;

    complex *R, *R2, *R_Com;

    MATRIX *fft_lhs, *fft_rhs;

    double real_z_a, imag_z_a;
    double real_z_b_con, imag_z_b_con;
    double real_a_b, imag_a_b;

    double real_z_c, imag_z_c;
    double real_z_d_con, imag_z_d_con;
    double real_c_d, imag_c_d;

    N= signal_length * 2;
	size_window = (window_length+1)/2;

    R =	(complex *) mxCalloc (window_r2, sizeof(complex));
    if (R==NULL)
    {
        mexPrintf( "\nMemory allocation	failed\n"	);
    }

    R2 =	(complex *) mxCalloc (window_r2, sizeof(complex));
    if (R2==NULL)
    {
        mexPrintf( "\nMemory allocation	failed\n"	);
    }

    R_Com =	(complex *) mxCalloc (window_r2, sizeof(complex));
    if (R_Com==NULL)
    {
        mexPrintf( "\nMemory allocation	failed\n"	);
    }

    fft_rhs = MXCREATEFULL( 1, window_r2, COMPLEX );

    data1 = (double *) mxMalloc(window_r2 * sizeof(double));
    data2 = (double *) mxMalloc(window_r2 * sizeof(double));

    data3 = (double *) mxMalloc(window_r2 * sizeof(double));
    data4 = (double *) mxMalloc(window_r2 * sizeof(double));

    z_real = (double *) mxMalloc(signal_length*2 * sizeof(double));
    z_imag = (double *) mxMalloc(signal_length*2 * sizeof(double));


    //---------------------------------//
    // Calling Analytic Signal Function//
    analytic_signal( signal_real, z_real, z_imag, signal_length);

	///////////////////////
    // Calculating Instantaneous Auto-correlation Function (IAF in Time-Lag Domain) //
    for (i = 0; i <=num_slices-1; (i=i+2))
    {
		t1=i*time_res;

		// Calculating IAF for each time-slice
        //R[j] =	z[t+j].z*[t-j],  with corrected indices:
        //
        for (j = 0; j <= size_window-1; j++)
        {
            f1=t1+j;
            f2=t1-j;

            // Computing Modulus for indexes
			a = Modulus_Index(f1, N);
			b = Modulus_Index(f2, N);

            real_z_a = z_real[a];
            imag_z_a = z_imag[a];

            real_z_b_con = z_real[b];
            imag_z_b_con = -z_imag[b];

            real_a_b = real_z_a*real_z_b_con - imag_z_a*imag_z_b_con;
            imag_a_b = real_z_a*imag_z_b_con + imag_z_a*real_z_b_con;

            R[j].re = G_real[0][j]*(real_a_b);
            R[j].im = G_real[0][j]*(imag_a_b);

            for (p = 1; p <= size_window-1; p++)
            {
                g1=p+f1;
                g2=p+f2;

                // Computing Modulus for indexes
                a = Modulus_Index (g1, N);
                b = Modulus_Index (g2, N);

                g3=-p+f1;
                g4=-p+f2;

                c = Modulus_Index (g3, N);
                d = Modulus_Index (g4, N);

                real_z_a = z_real[a];
                imag_z_a = z_imag[a];

                real_z_b_con = z_real[b];
                imag_z_b_con = -z_imag[b];

                real_z_c = z_real[c];
                imag_z_c = z_imag[c];

                real_z_d_con = z_real[d];
                imag_z_d_con = -z_imag[d];

                real_a_b = real_z_a*real_z_b_con - imag_z_a*imag_z_b_con;
                imag_a_b = real_z_a*imag_z_b_con + imag_z_a*real_z_b_con;

                real_c_d = real_z_c*real_z_d_con - imag_z_c*imag_z_d_con;
                imag_c_d = real_z_c*imag_z_d_con + imag_z_c*real_z_d_con;

                R[j].re = R[j].re + (G_real[p][j]*(real_a_b + real_c_d));
                R[j].im = R[j].im + (G_real[p][j]*(imag_a_b + imag_c_d));

            }

            if (j>0)
            {
                R[window_r2-j].re = R[j].re;
                R[window_r2-j].im = -R[j].im;
            }
        }

        t1=(i+1)*time_res;


        for (j = 0; j <= size_window-1; j++)
        {
            f1=t1+j;
            f2=t1-j;

			// Computing Modulus for indexes
            a = Modulus_Index (f1, N);
            b = Modulus_Index (f2, N);

            real_z_a = z_real[a];
            imag_z_a = z_imag[a];

            real_z_b_con = z_real[b];
            imag_z_b_con = -z_imag[b];

            real_a_b = real_z_a*real_z_b_con - imag_z_a*imag_z_b_con;
            imag_a_b = real_z_a*imag_z_b_con + imag_z_a*real_z_b_con;

            R2[j].re = G_real[0][j]*(real_a_b);
            R2[j].im = G_real[0][j]*(imag_a_b);

            for (p = 1; p <= size_window-1; p++)
            {
                g1=p+f1;
                g2=p+f2;

                // Computing Modulus for indexes
                a = Modulus_Index (g1, N);
                b = Modulus_Index (g2, N);

                g3=-p+f1;
                g4=-p+f2;

                // Computing Modulus for indexes
                c = Modulus_Index (g3, N);
                d = Modulus_Index (g4, N);

                real_z_a = z_real[a];
                imag_z_a = z_imag[a];

                real_z_b_con = z_real[b];
                imag_z_b_con = -z_imag[b];

                real_z_c = z_real[c];
                imag_z_c = z_imag[c];

                real_z_d_con = z_real[d];
                imag_z_d_con = -z_imag[d];

                real_a_b = real_z_a*real_z_b_con - imag_z_a*imag_z_b_con;
                imag_a_b = real_z_a*imag_z_b_con + imag_z_a*real_z_b_con;

                real_c_d = real_z_c*real_z_d_con - imag_z_c*imag_z_d_con;
                imag_c_d = real_z_c*imag_z_d_con + imag_z_c*real_z_d_con;

                R2[j].re = R2[j].re + (G_real[p][j]*(real_a_b + real_c_d));
                R2[j].im = R2[j].im + (G_real[p][j]*(imag_a_b + imag_c_d));

            }

            if (j>0)
            {
                R2[window_r2-j].re = R2[j].re;
                R2[window_r2-j].im = -R2[j].im;
			}
        }

		//---- Adding two time slices, by adding their real and imaginary parts
        for (ind = 0; ind <= window_r2-1; ind++)
        {
            R_Com[ind].re = R[ind].re - R2[ind].im;
            R_Com[ind].im = R[ind].im + R2[ind].re;

            data1[ind] = R_Com[ind].re;
            data2[ind] = R_Com[ind].im;
        }

	 	 // Take fourier transform in Lag domain to form Time-Frequency Distribution
         mxSetPr(fft_rhs, data1);
         mxSetPi(fft_rhs, data2);
         mexCallMATLAB(1, &fft_lhs, 1, &fft_rhs, "fft" );
         data3 = mxGetPr(fft_lhs);
         data4 = mxGetPi(fft_lhs);

         for (ind = 0; ind <= window_r2-1; ind++)
         {
             if (i<=num_slices-2)
             {
                result_real[ind+window_r2*i]=data3[ind];
                result_real[ind+window_r2*(i+1)]=data4[ind];

             }
             else
             {
                 result_real[ind+window_r2*i]=data3[ind];
             }
         }
    }

	mxDestroyArray(fft_rhs);

	return 0;
}

// Function for Calculating Modulus of Indexes
int Modulus_Index(int ind, int N)
{
    int out_ind;

    out_ind = ind % N;

    if (out_ind<=0)
    {
        out_ind=out_ind+N;
    }

    out_ind = out_ind-1;

    return out_ind;
}

// Function for Calculating Analytic Signal
void analytic_signal( double *signal_real,
        double *z_real_out, double *z_imag_out, int signal_length){

    int i;
    double *z_real, *z_imag;
    double *z_real_extend, *z_imag_extend;
    double *z_ifft_real, *z_ifft_imag;
    double *z_t_real, *z_t_imag;

    MATRIX *fft_lhs_data, *rhs_2_data, *fft_rhs_data, *lhs_2_data;

    rhs_2_data = MXCREATEFULL( 1, signal_length*2, COMPLEX );
    fft_rhs_data = MXCREATEFULL( 1, signal_length*2, COMPLEX );

    z_real = (double *) mxMalloc(signal_length*2 * sizeof(double));
    z_imag = (double *) mxMalloc(signal_length*2 * sizeof(double));
    z_real_extend = (double *) mxMalloc(signal_length*2 * sizeof(double));
    z_imag_extend = (double *) mxMalloc(signal_length*2 * sizeof(double));
    z_t_real = (double *) mxMalloc(signal_length*2 * sizeof(double));
    z_t_imag = (double *) mxMalloc(signal_length*2 * sizeof(double));
    z_ifft_real = (double *) mxMalloc(signal_length*2 * sizeof(double));
    z_ifft_imag = (double *) mxMalloc(signal_length*2 * sizeof(double));

    // Extending the length of input signal of length N to 2*N by zero padding
    for (i=0;i<signal_length*2;i++) {
        if (i<signal_length) {
            z_real_extend[i] = signal_real[i];
            z_imag_extend[i] = 0;
        }
        else {
            z_real_extend[i] = 0;
            z_imag_extend[i] = 0;
        }
    }

    mxSetPr(fft_rhs_data, z_real_extend);
    mxSetPi(fft_rhs_data, z_imag_extend);
    // Calling Matlab function fft for computing fourier transform of signal
    mexCallMATLAB(1, &fft_lhs_data, 1, &fft_rhs_data, "fft" );
    z_ifft_real = mxGetPr(fft_lhs_data);
    z_ifft_imag = mxGetPi(fft_lhs_data);

    // Multiplying the first half of the entexded input signal by 2
    for (i=1;i<signal_length;i++) {
        z_ifft_real[i] = 2*z_ifft_real[i];
        z_ifft_imag[i] = 2*z_ifft_imag[i];
    }
    // Replacing the second half of the entexded input signal by 0
    for (i=signal_length+1;i<2*signal_length;i++) {
        z_ifft_real[i] = 0;
        z_ifft_imag[i] = 0;
    }

    for (i=0;i<2*signal_length;i++) {
        z_t_real[i] = z_ifft_real[i];
        z_t_imag[i] = z_ifft_imag[i];
    }
    mxSetPr(rhs_2_data, z_t_real);
    mxSetPi(rhs_2_data, z_t_imag);
    // Calling Matlab function ifft for computing inverse fourier transform
    mexCallMATLAB(1, &lhs_2_data, 1, &rhs_2_data, "ifft" );
    z_real = mxGetPr(lhs_2_data);
    z_imag = mxGetPi(lhs_2_data);

    for (i=0;i<2*signal_length;i++) {
        z_real_out[i] = z_real[i];
        z_imag_out[i] = z_imag[i];
    }
}
