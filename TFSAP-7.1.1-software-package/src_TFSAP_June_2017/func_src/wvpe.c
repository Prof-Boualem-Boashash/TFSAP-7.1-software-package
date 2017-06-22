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
* from	the peak of the	Wigner-Ville Distribution.
*************************************************/

#include <math.h>
#include <ctype.h>
#include <stdlib.h>

#include "mex.h"

#include "arithm.h"		     /*	routines for complex data */
#include "tlocal.h"		     /*	local function prototypes */

#include "wvpe.h"
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



int wvpe( double *signal_real, double *signal_imag, int	signal_length,
	  double *result, int nplts, int time_res,
	  int window_length, int window_r2,
	  int window_order)
{
  int i, j, k, signal_length_r2, sliding_time;
  complex *sig,	*Wig, Kernel;
  double sampling_freq,	posmax;
  int hlf, i1, i2;
  double max;

  /*
   * Calculate radix-2 value above signal_length
   */
  signal_length_r2 = 1;
  while	(signal_length_r2 < signal_length)
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


  /* Allocate complex array storage for	the input signal */
  sig =	(complex *) mxCalloc (signal_length, sizeof(complex));
  if (sig==NULL)
  {
    tfsaErr( "wvpe", "Memory allocation	failed");
    return 1;
  }

  Wig =	(complex *) mxCalloc (window_r2, sizeof(complex));
  if (Wig == NULL)
  {
    tfsaErr( "wvpe", "Memory allocation	failed");
    return 1;
  }

  

  /* Read in input signal into complex array */
  if (signal_imag == NULL)
    {
      for (i = 0; i < signal_length; i++)
	{
	  sig[i].re = signal_real[i];
	  sig[i].im = 0.0;
	}

      /* Now generate analytic signal of the real signal */
      default_sigana(sig, signal_length);
    }
  else
    for	(i=0; i	< signal_length; i++)
    {
      sig[i].re	= signal_real[i];
      sig[i].im	= signal_imag[i];
    }


  sliding_time = 0;
  for (i=0; i<nplts; i+=2)
    {  /*
	*/

#ifdef WAITBAR
      *n1 = (double)i/nplts;
      mexCallMATLAB( 0,	plhs, 1, prhs, "waitbar");
#endif


      hlf = window_length/2;
      if ((int)sliding_time - (int)hlf < 0)
	{
	  hlf =	sliding_time;
	}
      if (sliding_time+hlf >= signal_length)
	{
	  hlf =	signal_length -1 -sliding_time;
	}

      for (j=hlf+1; j <	window_r2 - hlf; j++)
	{
	  Wig[j].re = 0.0;
	  Wig[j].im = 0.0;
	}

      for (j=0;	j <= hlf ; j++)
	{
	  /* Indices for two Kernel terms that occur at	integer	lags */
	  i1 = sliding_time + j;
	  i2 = sliding_time - j;

	  Kernel.re = sig[i1].re * sig[i2].re +	sig[i1].im * sig[i2].im;
	  Kernel.im = sig[i2].re * sig[i1].im -	sig[i1].re * sig[i2].im;


	  /* Prepare Kernel for	FFT */
	  Wig[j].re = Kernel.re;
	  Wig[j].im = Kernel.im;
	  if (j	!= 0)
	    {
	      Wig[window_r2-j].re = Kernel.re;
	      Wig[window_r2-j].im =-Kernel.im;
	    }
	}

      if (i < nplts-1)
	{
	  sliding_time += time_res;
	  hlf =	window_length/2;
	  if ((int)sliding_time	- (int)hlf < 0)
	    {
	      hlf = sliding_time;
	    }
	if (sliding_time+hlf >=	signal_length)
	  {
	    hlf	= signal_length	-1 -sliding_time;
	  }
	  for (j=0; j <= hlf ; j++)
	  {
	    /* Indices for two Kernel terms that occur at integer lags */
	    i1 = sliding_time +	j;
	    i2 = sliding_time -	j;

	    Kernel.re =	sig[i1].re * sig[i2].re	+ sig[i1].im * sig[i2].im;
	    Kernel.im =	sig[i2].re * sig[i1].im	- sig[i1].re * sig[i2].im;

	    /* Prepare Kernel for FFT */
	    Wig[j].re +=-Kernel.im;
	    Wig[j].im += Kernel.re;
	    if (j != 0)
	      {
		Wig[window_r2-j].re += Kernel.im;
		Wig[window_r2-j].im += Kernel.re;
	      }
	  }

	}


      FFT (Wig,	window_order, 1);


      sampling_freq = 1.;

      /* Find out the
       * frequency value at which WIG is Maximum.
       */
      max = Wig[0].re;
      posmax = 0.;
      for (k=1;	k<window_r2; k++)
	if (Wig[k].re >	max)
	  {
	    max	= Wig[k].re;
	    posmax = (double)k;
	  }
      *(result + i) = sampling_freq * posmax / (2.*(double) window_r2);

      if (i < nplts-1)
	{
	  /* Find out the
	   * frequency value at	which WIG is Maximum.
	   */
	  max =	Wig[0].im;
	  posmax = 0.;
	  for (k=1; k<window_r2; k++)
	    if (Wig[k].im > max)
	  {
	    max	= Wig[k].im;
	    posmax = (double)k;
	  }
	  *(result + i+1) = sampling_freq * posmax / (2.*(double) window_r2);
	}


      sliding_time+= time_res;
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








