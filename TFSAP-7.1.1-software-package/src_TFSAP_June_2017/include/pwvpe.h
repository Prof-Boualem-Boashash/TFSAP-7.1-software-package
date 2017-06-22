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
* Lastly the following reference is useful for understanding the basics of Instantaneous
* Frequency estimation:
* [4] B. Boashash, "Estimating and interpreting the instantaneous frequency of
* a signal—part 2: algorithms and applications", Proc. IEEE 80 (4) (1992) 540-568.
*
*
* Description:
 * Instantaneous Frequency Estimation by Peak of 
 * Polynomial Wigner-Ville Distribution
 ********************************************************************/

#ifndef ARITHM_H
#include "arithm.h"
#endif

/*
 * int wvd( double *signal_real, double *signal_imag, int
 *          signal_length, double *result, int num_slices,
 *          int time_res, int window_length, int
 *          window_r2, int window_order)
 *
 * Description:
 * 
 *     Computes the wvd of the input signal.  An analytic signal
 *     generator is called if the input signal is real.  The supplied
 *     length of the analysis window defines whether the "true" wvd or
 *     a pseudo (windowed) wvd is computed.  This version takes
 *     advantage of the reality of the wvd to increase calculation
 *     speed.
 *
 * Parameters:
 *
 *     double *signal_real
 *     double *signal_imag
 *
 *         Real and imaginary components of the input signals.  If
 *         signal_imag is passed as NULL, then the real input signal
 *         is made analytic within this function.  If signal_imag !=
 *         NULL, then it is assumed that the user is passing an
 *         analytic signal, and no futher pre-processing is done.
 *
 *     int signal_length
 *
 *         Length of signal_real and _imag arrays.
 *
 *     double *result
 *
 *         Should point to allocated memory to store the computed wvd.
 *         These result arrays are returned holding window_r2 *
 *         num_slices doubles each; num_slices of window_r2
 *         consecutive frequency values (beginning at 0 sec and 0 Hz
 *         samples).
 *
 *     int num_slices
 *
 *         The number of time instants at which the wvd will be
 *         computed.  This is closely related to time_res.  See below
 *         for more details and how to calculate it.
 *
 *     int time_res
 *
 *         The time-spacing between successive frequency slices of the
 *         wvd.  This is closely related to num_slices.  See below for
 *         more details.
 *
 *     int window_length
 *
 *         The time length of the window to be applied to the signal
 *         in generating the wvd.  Window_length must be odd to ensure
 *         that the wvd is not calculated between integer time
 *         samples.  If window_length is longer than or equal to the
 *         signal length, then the (normal) wvd is generated.  If the
 *         window_length is shorter, then the psuedo-wvd is generated
 *         (by definition).
 *
 *         NOTE: this function does not check for even window_length.
 *         It is the callers responsibility to comply.
 *
 *     int window_r2
 *
 *         The smallest power of two equal to or larger than
 *         window_length.  See below for rationale and calculation.
 *         Extra zero-padding at the fft stage of the analysis may be
 *         specified using a larger than normal value here.
 * 
 *     int window_order
 *
 *         Log2(window_r2).  See below for rationale and calculation.
 *
 * Comments:
 *
 *     Several of the parameters are redundant; specifically,
 *     num_slices can be determined from time_res, and window_r2,
 *     window_order can be determined from window_length.  However,
 *     since the calling program must work out the values of these
 *     parameters in order to correctly allocate memory for the
 *     result, they can be passed here to eliminate computing these
 *     values a second time.  To calculate the required values use the
 *     following code.
 *
 *     >     // begin calculate paramters required by wvd //
 *     >       signal_length = ...;
 *     >       time_res = ...;
 *     >       window_length = ...;
 *     >     
 *     >       num_slices = (int) ceil((double)signal_length/time_res);
 *     >     
 *     >       // calculate radix-2 value equal or above window_length //
 *     >       window_order = 0;
 *     >       window_r2 = 1;
 *     >       while( window_r2 < window_length) {
 *     >         window_order++;
 *     >         window_r2 <<= 1;
 *     >       }
 *     >     // end //
 *
 *     Memory for the results array must be allocated by the calling
 *     program using
 *
 *     >     // begin memory allocation //
 *     >       result = (double *)malloc( sizeof(double) * window_r2 * nplts);
 *     >     // end //
 *
 * Returns:
 *
 *     0       no errors
 *     1       memory allocation failed
 *
 */


void pwvpe( double *signal_real, double *signal_imag, int signal_length, 
          double *result, int nplts, int time_res, 
          int window_length, int window_r2, 
          int window_order, int deg);


