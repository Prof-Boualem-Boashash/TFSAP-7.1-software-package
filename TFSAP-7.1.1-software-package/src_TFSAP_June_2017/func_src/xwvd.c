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
* Program for generating the Cross Wigner-Ville and Pseudo Cross Wigner-Ville
* Distributions.  This	is based on wvd.
*************************************************/

#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mex.h"


#include "arithm.h"		     /* routines for complex data */
#include "tlocal.h"		     /* local function prototypes */

#include "analyt.h"
#include "fft.h"
#include "xwvd.h"
#include "tfsa_c.h"

/* Global variables for waitbar */
#ifdef WAITBAR
MATRIX *prhs[2];
MATRIX *p1, *p2;
MATRIX *plhs[1];
double *n1;
double waitbar_hndl;
#endif

int xwvd( double *signal1_real, double *signal1_imag, double
         *signal2_real, double *signal2_imag, int signal_length, double
	 *result_real, double *result_imag, int num_slices, int
	 time_res, int window_length, int window_r2, int
	 window_order)
{
  int i, j;
  complex *sig1, *sig2, *wvd;
  int w_centre;
  int hlf;
  int i1, i2;

  /* read in input signal into complex array */
  sig1 = (complex *) mxCalloc (signal_length, sizeof(complex));
  if (sig1 == NULL) {
    tfsaErr( "xwvd", "Memory allocation failed");
    return 1;
  }
  sig2 = (complex *) mxCalloc (signal_length, sizeof(complex));
  if (sig2 == NULL) {
    tfsaErr( "xwvd", "Memory allocation failed");
/*    mxFree((void *)sig1); */
    return 1;
  }

  if (signal1_imag == NULL) {
    for (i = 0; i < signal_length; i++) {
      sig1[i].re = signal1_real[i];
      sig1[i].im = 0.0;
    }
    default_sigana(sig1,signal_length);
  }
  else
    for (i=0; i < signal_length; i++) {
      sig1[i].re = signal1_real[i];
      sig1[i].im = signal1_imag[i];
    }

  if (signal2_imag == NULL) {
    for (i = 0; i < signal_length; i++) {
      sig2[i].re = signal2_real[i];
      sig2[i].im = 0.0;
    }
    default_sigana(sig2,signal_length);
  }
  else
    for (i=0; i < signal_length; i++) {
      sig2[i].re = signal2_real[i];
      sig2[i].im = signal2_imag[i];
    }

  /* allocate storage for slice of wvd */
  wvd = (complex *) mxCalloc (window_r2, sizeof(complex));
  if (wvd == NULL) {
    tfsaErr( "xwvd", "Memory allocation failed");
    /* mxFree((void *)sig1);
    mxFree((void *)sig2); */
    return 1;
  }

  /* starting at 0 is not real clever, since the first slice will
     apply a window of 1 long. For want of a better procedure, we
     will start at the first time_res step.  This is why we have the
     minus one at the end of the line to compute the number of slices.
     See function header comments for more info. */
  w_centre = 0;


#ifdef WAITBAR
  p1 = mxCreateFull(1, 1, REAL);
  p2 = mxCreateString("Please wait, computing distribution...");
  n1 = mxGetPr(p1);
  *n1 = 0.0;
  prhs[0] = p1;
  prhs[1] = p2;
  mexCallMATLAB( 1, plhs, 2, prhs, "waitbar");
  waitbar_hndl = *mxGetPr(plhs[0]);
#endif


  for (i = 0; i < num_slices; i++) {

#ifdef WAITBAR
    *n1 = (double)i/num_slices;
    mexCallMATLAB( 0, plhs, 1, prhs, "waitbar");
#endif

    /* now adjust effective window length to ensure that window is
       contained within signal */

    /* note that window length across slices (for 2d kernel) will be
       the same as the window length along kernel */

    hlf = window_length/2;
    if ((int)w_centre-(int)hlf < 0) {
      hlf = w_centre;
      /* effective_wl = w_centre * 2 + 1; */
    }
    if (w_centre+hlf>=signal_length) {
      hlf = signal_length-1-w_centre;
      /* effective_wl = (signal_length-1-w_centre) * 2 + 1; */
    }

    /* initialize wvd array (only those points not explicitly set
       below */
    for (j=hlf+1; j < window_r2-hlf; j++) {
      wvd[j].re = 0.0;
      wvd[j].im = 0.0;
    }

    for (j=0; j <= hlf ; j++) {
      i1 = w_centre + j;
      i2 = w_centre - j;
      wvd[j].re = sig1[i1].re * sig2[i2].re + sig1[i1].im * sig2[i2].im;
      wvd[j].im = sig2[i2].re * sig1[i1].im - sig1[i1].re * sig2[i2].im;
      if (j != 0) {
        wvd[window_r2-j].re = sig1[i2].re * sig2[i1].re + sig1[i2].im * sig2[i1].im;
	wvd[window_r2-j].im = sig2[i1].re * sig1[i2].im - sig1[i2].re * sig2[i1].im;
      }
    }

    FFT(wvd, window_order, 1);

    if (result_imag == NULL)
      for(j=0; j < window_r2; j++)
	*(result_real+i*window_r2+j) = sqrt( wvd[j].re * wvd[j].re + wvd[j].im * wvd[j].im);
    else
      for(j=0; j < window_r2; j++) {
	*(result_real+i*window_r2+j) = wvd[j].re;
	*(result_imag+i*window_r2+j) = wvd[j].im;
      }

    w_centre += time_res;
  }

#ifdef WAITBAR
  *n1 = (double)1;
  mexCallMATLAB( 0, plhs, 1, prhs, "waitbar");
  *n1 = waitbar_hndl;
  mexCallMATLAB( 0, plhs, 1, prhs, "close");
  mxFreeMatrix(p1);
  mxFreeMatrix(p2);
#endif

  return 0;
}

