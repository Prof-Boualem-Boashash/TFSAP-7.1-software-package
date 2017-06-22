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
* Routine for estimating the Instantaneous Frequency of a signal
* by adaptive recursive least square method
*************************************************/
#include <math.h>
#include <ctype.h>
#include <stdlib.h>

#include "mex.h"


#include "arithm.h"		     /*	routines for complex data types	*/
#include "tlocal.h"		     /*	local function prototypes */

#include "rls.h"
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


int rls(double *signal_real, double *signal_imag, int signal_length,
	double *result,	double alpha)
{
   int i, n;
   complex *sig, pre_coeff, error, Psig;
   double IFE, P;

   /* Allocate complex array storage for the input signal */
   sig = (complex *) mxCalloc (signal_length, sizeof(complex));
   if( sig == NULL ) {
      tfsaErr( "rms", "Internal	memory failure:	sig" );
      return 1;
   }


   /* Read in input signal into	complex	array */
   if (signal_imag == NULL)
     {
	for (i = 0; i <	signal_length; i++)
	  {
	     sig[i].re = signal_real[i];
	     sig[i].im = 0.0;
	  }

	/* Now generate	analytic signal	of the real signal */
	default_sigana(sig,signal_length);
     }
   else
     for (i=0; i < signal_length; i++)
     {
	sig[i].re = signal_real[i];
	sig[i].im = signal_imag[i];
     }


   /*
    * Initialize Waitbar
    */

#ifdef WAITBAR
   p1 =	mxCreateFull(1,	1, REAL);
   p2 =	mxCreateString("Please wait, estimating	frequency...");
   n1 =	mxGetPr(p1);
   *n1 = 0.0;
   prhs[0] = p1;
   prhs[1] = p2;
   mexCallMATLAB( 1, plhs, 2, prhs, "waitbar");
   waitbar_hndl	= *mxGetPr(plhs[0]);
#endif


   P = 1.0;
   pre_coeff.re	= sig[0].re;
   pre_coeff.im	= sig[0].im;
   for (n = 1; n < signal_length; n++)
     {
	error.re = pre_coeff.re*sig[n].re - pre_coeff.im*sig[n].im + sig[n+1].re;
	error.im = pre_coeff.re*sig[n].im + pre_coeff.im*sig[n].re + sig[n+1].im;

	P = (alpha/P) +	sig[n].re*sig[n].re + sig[n].im*sig[n].im;
	P = 1.0	/ P;

	Psig.re	= P * sig[n].re;
	Psig.im	= P * sig[n].im;

	pre_coeff.re = pre_coeff.re -
	  (Psig.re * error.re +	Psig.im	* error.im);
	pre_coeff.im = pre_coeff.im -
	  (Psig.re * error.im -	Psig.im	* error.re);

	IFE = atan(pre_coeff.im/pre_coeff.re);
	if (IFE	< 0) IFE += TFSA_PI;
	*(result + n) =	IFE / (2.0*TFSA_PI);

#ifdef WAITBAR
	*n1 = (double)n/signal_length;
	mexCallMATLAB( 0, plhs,	1, prhs, "waitbar");
#endif

     }


#ifdef WAITBAR
   *n1 = (double)1;
   mexCallMATLAB( 0, plhs, 1, prhs, "waitbar");
   *n1 = waitbar_hndl;
   mexCallMATLAB( 0, plhs, 1, prhs, "close");
   mxFreeMatrix(p1);
   mxFreeMatrix(p2);
#endif


   mxFree((void	*)sig);

   return 0;
}
