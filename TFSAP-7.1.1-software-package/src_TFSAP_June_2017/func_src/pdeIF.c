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
* [2] B. Boashash, Samir Ouelha, Efficient Software Matlab Package for the calculation of Time-Frequency Distributions related
* Time-Scale methods and the extraction of signal characteristics, SoftwareX, 2017.
*
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
* by General Phase Difference (for orders 1,2,4)  method and the
* weighted phase difference  Kay  estimator using a second order
* Kay window
*************************************************/
#include <math.h>
#include <ctype.h>
#include <stdlib.h>

#include "mex.h"

#include "arithm.h"		     /*	routines for complex data */
#include "tlocal.h"		     /*	local function prototypes */

#include "pdeIF.h"
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


int pdeIF(double *signal_real, double *signal_imag, int signal_length,
	double *result,	int window_length, int estimator_order,
	int Kay)
{
  double *IFE, temp, temp1, temp2, smooth;
  int i, j, signal_length_r2;
  complex *sig,	X;
  double *unwrapped;
  double Phase,	phase, phase1;
  double sum;
  void unwrapping ( complex *sig, int signal_length, double *temp);



  /*
   * Calculate radix-2 value above signal_length
   */
   signal_length_r2 = 1;
   while (signal_length_r2 < signal_length)
     signal_length_r2 <<= 1;

  /*
   * Initialize	Waitbar
   */

#ifdef WAITBAR
  p1 = mxCreateFull(1, 1, REAL);
  p2 = mxCreateString("Please wait, estimating frequency...");
  n1 = mxGetPr(p1);
  *n1 =	0.0;
  prhs[0] = p1;
  prhs[1] = p2;
  mexCallMATLAB( 1, plhs, 2, prhs, "waitbar");
  waitbar_hndl = *mxGetPr(plhs[0]);
#endif

#ifdef MATLAB_MEX_FILE
  IFE =	(double	*) mxCalloc (signal_length -window_length +1, sizeof(double));
   if (IFE == NULL)
     {
	tfsaErr("pde", "Memory allocation failed" );
	return 1;
     }
#else
   IFE = (double *) calloc (signal_length -window_length +1, sizeof(double));
   if (IFE == NULL)
     {
	tfsaErr("pde", "Memory allocation failed" );
	return 1;
     }

#endif

  /* Allocate complex array storage for	the input signal */
#ifdef MATLAB_MEX_FILE
  sig =	(complex *) mxCalloc (signal_length, sizeof(complex));
  if (sig==NULL) {
     tfsaErr("pde", "Memory allocation failed" );
     return 1;
  }
#else
  sig =	(complex *) calloc (signal_length, sizeof(complex));
  if (sig==NULL)
  {
     tfsaErr("pde", "Memory allocation failed" );
     return 1;
  }
#endif

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

  if (Kay)
  {
    /* Now smooth the phase of the analytic signal using a
     * quadratic Kay window
     */
    for	(i=0; i<signal_length-window_length; i++)
    { /* sliding window	a long the signal */

#ifdef WAITBAR
      *n1 = (double)i/(signal_length-window_length);
      mexCallMATLAB( 0,	plhs, 1, prhs, "waitbar");
#endif

       /* Now compute the local	smoothed estimate inside the
	* moving quadratic window each time the	window moves
	*/
       sum = 0.0;
       for (j=0; j<=window_length-2; j++)
       {
	  X.re = sig[i+j].re * sig[i+j+1].re + sig[i+j].im * sig[i+j+1].im;
	  X.im = sig[i+j].re * sig[i+j+1].im - sig[i+j].im * sig[i+j+1].re;
	  Phase	= atan(X.im/X.re);

	  /*Phase=atan(sig[i+j+1].im/sig[i+j+1].re)-atan(sig[i+j].im/sig[i+j].re);
*/
	  if (Phase < (double)0) Phase += TFSA_PI;

	  /* compute smoothing coefficient */
	  temp = (double) window_length;
	  temp1	=temp/2.;
	  temp2	= ((double)(j) - (temp1-1.))/temp1;
	  temp2	= temp2	* temp2;
	  smooth = (1.5*temp/(temp*temp-1)) * (1.-temp2);

	  sum += (smooth * (double)Phase);
       }

       /* Now place each local estimated frequency at its
	* window center
	*/
       *(result	+ window_length/2 + i) = sum/(2.*TFSA_PI);
    }
  }

  else	/* Generalized Finite Difference Estimator */
  {

    if ((estimator_order == 4)||(estimator_order == 6))
      {
#ifdef MATLAB_MEX_FILE
  unwrapped = (double *) mxCalloc (signal_length, sizeof(double));
  if (unwrapped	== NULL)
  {
     tfsaErr("pde", "Memory allocation failed" );
     return 1;
  }
#else
  unwrapped = (double *) calloc	(signal_length,	sizeof(double));
  if (unwrapped	== NULL)
  {
     tfsaErr("pde", "Memory allocation failed" );
     return 1;
  }
#endif
	unwrapping(sig,signal_length,unwrapped);
      }


    for	(i=estimator_order/2; i<signal_length-(estimator_order/2); i++)
    {

#ifdef WAITBAR
      *n1 = (double)i/(signal_length-(estimator_order));
      mexCallMATLAB( 0,	plhs, 1, prhs, "waitbar");
#endif

       switch (estimator_order){
      case 1: /* 1st order estimator (FFD) */
	if (i<signal_length -1)
	   Phase = atan(sig[i+1].im/sig[i+1].re)-atan(sig[i].im/sig[i].re);
	   if (Phase < (double)0) Phase	+= TFSA_PI;
	break;

      case 2: /* 2nd order estimator (CFD) */
	X.re = sig[i].re * sig[i+1].re + sig[i].im * sig[i+1].im;
	X.im = sig[i].re * sig[i+1].im - sig[i].im * sig[i+1].re;
	phase =	atan(X.im/X.re);
	if (phase < (double)0) phase +=	TFSA_PI;
	Phase =	phase;
	X.re = sig[i].re * sig[i-1].re + sig[i].im * sig[i-1].im;
	X.im = sig[i].im * sig[i-1].re - sig[i-1].im * sig[i].re;
	phase =	atan(X.im/X.re);
	if (phase < (double)0) phase +=	TFSA_PI;
	Phase += phase;
	Phase *= 0.5;
	break;

      case 4: /* 4th order estimator (with phase unwrapping) */
	phase =	(double)unwrapped[i-2];
	Phase =	phase;
	phase =	(double)unwrapped[i+2];
	Phase += -phase;
	Phase *= (double)(1./12.);

	X.re = sig[i].re * sig[i+1].re + sig[i].im * sig[i+1].im;
	X.im = sig[i].re * sig[i+1].im - sig[i].im * sig[i+1].re;
	phase1 = atan(X.im/X.re);
	if (phase1 < (double)0)	phase1 += TFSA_PI;
	X.re = sig[i].re * sig[i-1].re + sig[i].im * sig[i-1].im;
	X.im = sig[i].im * sig[i-1].re - sig[i-1].im * sig[i].re;
	phase =	atan(X.im/X.re);
	if (phase < (double)0) phase +=	TFSA_PI;
	phase += phase1;
	phase *= (double)(2./3.);
	Phase += phase;

	break;

      case 6: /* 6th order estimator (with phase unwrapping) */
	phase =	unwrapped[i+3];
	Phase =	phase;
	phase =	unwrapped[i-3];
	Phase += -phase;
	Phase *= (double)(1./60.);

	phase =	unwrapped[i-2];
	phase1 = unwrapped[i+2];
	phase += -phase1;
	phase *= (double)(3./20.);
	Phase += phase;

	X.re = sig[i].re * sig[i+1].re + sig[i].im * sig[i+1].im;
	X.im = sig[i].re * sig[i+1].im - sig[i].im * sig[i+1].re;
	phase1 = atan(X.im/X.re);
	if (phase1 < (double)0)	phase1 += TFSA_PI;
	X.re = sig[i].re * sig[i-1].re + sig[i].im * sig[i-1].im;
	X.im = sig[i].im * sig[i-1].re - sig[i-1].im * sig[i].re;
	phase =	atan(X.im/X.re);
	if (phase < (double)0) phase +=	TFSA_PI;
	phase += phase1;
	phase *= (double)0.75;
	Phase += phase;

	break;

	default:
	  tfsaErr( "pde", "Illegal estimator order..." );
	  return 1;
	  break;
       }

      *(result + i) = (double) Phase/(2.*TFSA_PI);
    }
  }

  if ((estimator_order == 4)||(estimator_order == 6))
    {


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


void unwrapping	( complex *sig,	int signal_length, double *temp)
/*
 * This	routine	unwrap the phase of the	signal sig
 */
{
  int i;
  double ren, rep, imn,	imp, g,	adj;


  for (i=0; i<signal_length; i++)
    {
      if (sig[i].re < 0)
	{
	  ren =	1;
	  rep =	0;
	}
      else if (sig[i].re > 0)
	{
	  ren =	0;
	  rep =	1;
	}
      if (sig[i].im < 0)
	{
	  imn =	1;
	  imp =	0;
	}
      else if (sig[i].im > 0)
	{
	  imn =	0;
	  imp =	1;
	}
      if (sig[i].re != 0)
	temp[i]	= atan(sig[i].im/sig[i].re) + ren*TFSA_PI + rep*imn*2.*TFSA_PI;
      else
	temp[i]	= imp*TFSA_PI/2. + imn*TFSA_PI*1.5;
      if (i > 0)
	{
	  adj =	temp[i-1] - g*2.*TFSA_PI;
	  if ((temp[i] < adj-TFSA_PI/2.) && (adj > TFSA_PI/2.))
	    g++;
	  if ((adj < TFSA_PI/2.) && (temp[i] > adj+1.5*TFSA_PI))
	      g--;
	  temp[i] += TFSA_PI*g*2.;
	}
    }
}
