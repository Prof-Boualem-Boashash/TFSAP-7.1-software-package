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
* The following references 2 should be cited whenever this script is used:
* [1] B. Boashash, Samir Ouelha, Designing time-frequency  and time-scale
* features for efficient classification of non stationary signals: a tutorial
* review with performance comparison, Digital Signal Processing, In Press.
* [2] B. Boashash, Samir Ouelha, Efficient Software Matlab Package for the calculation 
* of Time-Frequency Distributions related Time-Scale methods and the extraction of 
* signal characteristics, SoftwareX, In Press.
* In addition, the following 3rd reference is relevant for use of the executable
* [3] B. Boashash (ed.), Time-Frequency Signal Analysis and Processing, 2nd Edition 
* (London: Elsevier / Academic Press, December 2015); ISBN 978-0-12-398499-9. Book website:
* http://booksite.elsevier.com/9780123984999/
*
*
* Description: * Synthesize signal from a given TFD
 * 
 ********************************************************************
 */

#ifndef ARITHM_H
#include "arithm.h"
#endif


#define IDFT    1
#define OLA     2
#define MSTFT   3
#define MSPEC   4
#define MWVD    5

#define TOL_DEFAULT   1000

void mstft(double *y_real, double *y_imag, int M, double *synth_signal_re,
	   double *synth_signal_im, int signal_length, double *window,
	   int window_length, int hlf_win);


int inverseFFT(double *tfd_re, double *tfd_im, int tfd_M, int tfd_N, int M,
		double *y_real, double *y_imag);


void generateRandomSignal(double *synth_signal_re, double *synth_signal_im, 
			  int signal_length);


int do_svd( double *A_re, double *A_im, int a_M, double *u_re, double *u_im,
	    double *s_re, double *s_im );


int interpolate( double *synth_signal_re, double *synth_signal_im, int signal_length, 
		 complex *synth_sig_even );


int reconstructPhase( double *signal_re, double *signal_im, int signal_length, 
		      double *synth_signal_re, double *synth_signal_im );


int synthesize( int analysis_type, double *synth_signal_re, double *synth_signal_im,
		int signal_length, int window_length, int window_type, double tol, 
		double *signal_real, double *signal_imag, double *tfd_re, 
		double *tfd_im, int tfd_M, int tfd_N, double *rankone_err );







