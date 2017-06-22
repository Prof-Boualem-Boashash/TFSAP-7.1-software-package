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
* from	the peak of the	Time-Varying spectrogram.
*************************************************/

#include <math.h>
#include <ctype.h>
#include <stdlib.h>

#include "mex.h"


#include "arithm.h"		     /*	routines for complex data */
#include "tlocal.h"		     /*	local function prototypes */

#include "sfpe.h"
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



int sfpe( double  *signal_real, 
          double  *signal_imag, 
          int	  signal_length,
	  double  *result, 
          int     nplts, 
          int     time_res,
	  int     window_length, 
          int     window_r2,
	  int     window_order)
{
  int i, j, signal_length_r2, sliding_time;
  complex *sig,	*stft;
  double sampling_freq,	posmax;
  int hlf;
  double max;
  double mag;                                /* magnitude of the fft */


  /*
   * Calculate radix-2 value above signal_length
   */
  signal_length_r2 = 1;
  while	(signal_length_r2 < signal_length)
	signal_length_r2 <<= 1;


  /* Allocate complex array storage for	the input signal */
  sig =	(complex *) mxCalloc (signal_length, sizeof(complex));
  if (sig==NULL) {
    tfsaErr( "sfpe", "Memory allocation	failed");
    return 1;
  }

  stft = (complex *) mxCalloc (window_r2, sizeof(complex));
  if (stft == NULL)
  {
    tfsaErr( "sfpe", "Memory allocation	failed");
    return 1;
  }

  /* Read in input signal into complex array */
   if (signal_imag == NULL) {
      for (i = 0; i < signal_length; i++) {
         sig[i].re = signal_real[i];
         sig[i].im = 0.0;
      }

     /* Now generate analytic signal of the real signal */
     default_sigana(sig, signal_length);
   } else {
      for (i=0; i < signal_length; i++) {
         sig[i].re = signal_real[i];
         sig[i].im = signal_imag[i];
     }
   }
   /* Calculate	Spectogram and find out	its peaks */
   sliding_time	= time_res;
   for (i=0; i<nplts; i++) {

      /* Check	for window too long */
      for (j=0; j<window_r2; j++) {
	  stft[j].re =	0.;
          stft[j].im =	0.;
      }

      /* Now form Kernel for STFT */
      if (window_length > signal_length) {
         window_length	= signal_length;
      }

      hlf = window_length/2;
      if (sliding_time	< hlf) {
         hlf = sliding_time;
      }
      if (sliding_time	>= signal_length-hlf) {
         hlf = signal_length -sliding_time;
      }
      /* if (hlf == 0) return 1; */
      for (j=sliding_time-hlf; j<sliding_time+hlf; j++) {
         stft[j-sliding_time+hlf] = cmplx(sig[j].re, sig[j].im);
      } 
      /* Compute STFT for each	time slice and find out	the
       * frequency value at which STFT	is Maximum.
       */

      FFT( stft, window_order,	1);

      sampling_freq = 1.;
      max = stft[0].re*stft[0].re+stft[0].im*stft[0].im;
      for (j=1; j < window_r2; j++) {
        mag = stft[j].re*stft[j].re+stft[j].im*stft[j].im;
	if (mag > max) {
	   max = mag;
	   posmax = (double)j;
	}
      }
      *(result	+ i) = sampling_freq * posmax /	((double) window_r2);


      sliding_time += time_res;
   }


return 0;
}



void maxfft (complex *A, int M,	int N,
	 double	*valmax, double	*posn)
{
  int i,j;
  double coarse_max, coarse_ind;
  double magn,left_mag,	right_mag;
  double highind, botind, centind, valmid;
  double expn1,	expn;
  complex sum1,	sum, ctemp,ctemp1;
  complex *B;

  B = (complex *) mxCalloc (N, sizeof(complex));
  if (B	== NULL)
  {
    tfsaErr( "sfpe", "Memory allocation	failed");
    return;
  }

  for (i= 0 ; i < N; i++) { 
     B[i] = cmplx(A[i].re, A[i].im); 
  }

  FFT (A, M, 1);

  coarse_max = 0.;
  for (i=0; i <	N; i++) {
      magn = module(A[i]);
      if (magn > coarse_max) {
	   coarse_max =	magn;
	   coarse_ind =	(double) i;
      }
   }
   left_mag = module ( A [(int)coarse_ind - 1]);
   right_mag = module (	A [(int)coarse_ind + 1]	);
   *valmax = (double) coarse_max;
   highind = (double) coarse_ind;

   valmid = 0.;
   if (left_mag	> right_mag) { 
      botind = centind = (double)coarse_ind - 1.0;
   } else { 
      botind = centind = (double)coarse_ind + 1.;
   }

    for	(j = 1;	j <=10;	j++) {
	if (valmid > *valmax) {
	     botind = highind;
	     *valmax = valmid;
	     highind = centind;
	     centind = (highind+botind)/2.;
	 } else {
	     botind = centind;
	     centind = (highind + botind) /2.;
	 }

	 sum = cmplx(0., 0.);
	 sum1 =	cmplx(0., 0.);
	 for (i=0; i < N; i++) {
	      expn1 = 2.*TFSA_PI*(double)i*centind/(double)N;
	      expn = 2.*TFSA_PI*(double)i*coarse_ind/(double)N;
	      ctemp1 = cmplx(COS(expn1), -SIN(expn1));
	      ctemp  = cmplx(COS(expn),	-SIN(expn));
	      ctemp1 = multpl(ctemp1,B[i]);
	      ctemp  = multpl(ctemp, B[i]);
	      sum1  = add (sum1, ctemp1);
	      sum   = add (sum,	ctemp);
	  }
	  valmid = module ( sum1);
      }

      *valmax =	valmid;
      *posn = centind;

  }

