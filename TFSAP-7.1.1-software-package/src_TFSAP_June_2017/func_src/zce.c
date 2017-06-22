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
* using the Zero-Crossing method
*************************************************/

#include <math.h>
#include <ctype.h>
#include <stdlib.h>

#include "mex.h"

#include "arithm.h"		     /* routines for complex data types	*/
#include "tlocal.h"		     /* local function prototypes */

#include "zce.h"
#include "tfsa_c.h"

#ifdef WAITBAR
MATRIX *prhs[2];
MATRIX *p1, *p2;
MATRIX *plhs[1];
double *n1;
double waitbar_hndl;
#endif

#define TFSA_PI 3.14159265358979323846



int zce(double *signal_real, double *signal_imag, unsigned signal_length,
      double *result, unsigned window_length)
{
  unsigned n, m, ind, limit, half_window_length;
  double *sig;
  double Zero;
  int signum (double y);



   /* Allocate array storage for the input signal */
   sig = (double *) mxCalloc (signal_length, sizeof(double));
   if( sig == NULL ) {
      tfsaErr( "zce", "Memory allocation failed: sig" );
      return 1;
   }

  /* Read in input signal into an array */
  if (signal_imag == NULL)
    for (n = 0; n < signal_length; n++)
      sig[n] = signal_real[n];
  else
    for (n = 0; n < signal_length; n++)
      sig[n] = signal_real[n];


  /*
   * Initialize Waitbar
   */

#ifdef WAITBAR
  p1 = mxCreateFull(1, 1, REAL);
  p2 = mxCreateString("Please wait, estimating frequency...");
  n1 = mxGetPr(p1);
  *n1 = 0.0;
  prhs[0] = p1;
  prhs[1] = p2;
  mexCallMATLAB( 1, plhs, 2, prhs, "waitbar");
  waitbar_hndl = *mxGetPr(plhs[0]);
#endif


  half_window_length = window_length/2;
  for (n = 1 ; n < (signal_length-1); n++)
  {

/*     if (n%(signal_length/8)==0)*/
#ifdef WAITBAR
      *n1 = (double)n/(signal_length-2);
      mexCallMATLAB( 0, plhs, 1, prhs, "waitbar");
#endif

     Zero = 0.;
     if (n < half_window_length)
     {
         limit = 2*n;
	 for (m = 1; m <=limit; m++)
	    Zero = Zero+(double)abs( signum(sig[m]) - signum(sig[m-1]) );
	 Zero = Zero/(double)(4.*n);
     }
     else
     {
	 if (n < (signal_length - half_window_length ) )
	 {
	    limit = half_window_length +half_window_length;
	    for (m = 1; m <= limit; m++)
	    {
	      ind = n+m-half_window_length;
	      Zero = Zero+(double)abs( signum(sig[ind]) - signum(sig[ind-1]) );
	    }
	    Zero = Zero / (double)(4.*half_window_length );
	 }
	 else
	 {
	    limit = 2*(signal_length-n-1);
	    for (m = 1; m <= limit; m++)
	    {
	      ind = n+m-(signal_length-n-1);
	      Zero = Zero+(double)abs( signum(sig[ind]) - signum(sig[ind-1]) );
	    }
	    Zero = Zero / (double)(4. * (signal_length-n-1));
	 }
     }
     *(result+n) = (double)(Zero/2.);
  }

#ifdef WAITBAR
  *n1 = (double)1;
  mexCallMATLAB( 0, plhs, 1, prhs, "waitbar");
  *n1 = waitbar_hndl;
  mexCallMATLAB( 0, plhs, 1, prhs, "close");
  mxFreeMatrix(p1);
  mxFreeMatrix(p2);
#endif

  mxFree((void *)sig);

  return 0;
}


int signum (double y)
{
return(y<(double)0?-1:1);
}
