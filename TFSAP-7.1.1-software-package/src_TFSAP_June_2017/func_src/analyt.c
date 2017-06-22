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
*
* Description:
*
* Functions for the generation	of analytic signals from non-analytic signals
*
* For discrete case, Z[n] at n = Nyquist term (N/2) equals X[n],
* where X[.] = FT{ real signal x[.]) and
*       Z[.] = FT{ analytic signal z[.])
*
*
*
*************************************************/

#include <stdlib.h>
#define	huge __huge
#include "mex.h"

#include "arithm.h"		     /*	routines for complex data */
#include "tlocal.h"		     /*	local function prototypes */

#include "fft.h"
#include "analyt.h"
#include "tfsa_c.h"


void default_sigana( complex *X, unsigned signal_length)
{
  unsigned i;
  complex *Y;
  unsigned signal_order, signal_r2;

  signal_order = 0;
  signal_r2 = 1;
  while	(signal_r2 < signal_length) {
    signal_order++;
    signal_r2 <<= 1;
  }

  /* zero pad data array appropriately,	as required by FFT routine */
  Y = (complex *)mxCalloc( signal_r2, sizeof(complex));
  if (Y==NULL) {
    tfsaErr( "analyt", "Memory allocation failed: Y");
    return;
  }
  for (i = 0; i< signal_length;	i++) {
    Y[i].re = X[i].re;
    Y[i].im = X[i].im;
  }

  FFT(Y, signal_order, 1);
  /* zero negative frequencies */
  for (i=(signal_r2/2)+1; i<signal_r2;i++) {
    Y[i].re = 0.0;
    Y[i].im = 0.0;
  }
  /* multiply by two */
  for (i=1;   i	< signal_r2/2; i++) {
    Y[i].re += Y[i].re;
    Y[i].im += Y[i].im;
  }
  FFT(Y, signal_order, -1);

  /* discard zero-padded values.  Note that the	discarded values are
     not necessarily equivalent	to zero, and this truncation may be
     responsible for minor effects. */
  for (i = 0; i< signal_length;	i++) {
    X[i].re = Y[i].re;
    X[i].im = Y[i].im;
  }



}

