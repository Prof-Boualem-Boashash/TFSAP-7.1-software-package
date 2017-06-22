
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
* Description:
 *
 * xwvd
 *
 * Program for generating the Cross Wigner-Ville and Pseudo Cross Wigner-Ville
 * Distributions.  This is based on wvd.
 *
 *
 **/

/****
 *
 * int xwvd( double *signal1_real, double *signal1_imag, double
 * 	 *signal2_real, double *signal2_imag, unsigned signal_length, double
 * 	 *result_real, double *result_imag, unsigned num_slices, unsigned
 * 	 time_res, unsigned window_length, unsigned window_r2, unsigned
 * 	 window_order)
 *
 * Description:
 * 
 *     Computes the cross-wvd of the input signal.  An analytic signal
 *     generator is called if the input signals are real.  The supplied
 *     length of the analysis window defines whether the "true" wvd or
 *     a pseudo (windowed) (cross)-wvd is computed.
 *
 * Parameters:
 *
 *     double *signal1_real
 *     double *signal1_imag
 *     double *signal2_real
 *     double *signal2_imag
 *
 *         Real and imaginary components of the input signals.  If
 *         signal?_imag is passed as NULL, then the real input signal
 *         is made analytic within this function.  If signal?_imag !=
 *         NULL, then it is assumed that the user is passing an
 *         analytic signal, and no futher pre-processing is done.
 *
 *     unsigned signal_length
 *
 *         Length of signal?_real and _imag arrays.
 *
 *     double *result_real
 *     double *result_imag
 *
 *         Should point to allocated memory to store the computed wvd.
 *         These result arrays are returned holding window_r2 *
 *         num_slices doubles each; num_slices of window_r2
 *         consecutive frequency values (beginning at 0 sec and 0 Hz
 *         samples).
 *
 *     unsigned num_slices
 *
 *         The number of time instants at which the wvd will be
 *         computed.  This is closely related to time_res.  See below
 *         for more details and how to calculate it.
 *
 *     unsigned time_res
 *
 *         The time-spacing between successive frequency slices of the
 *         wvd.  This is closely related to num_slices.  See below for
 *         more details.
 *
 *     unsigned window_length
 *
 *         The time length of the window to be applied to the signal
 *         in generating the xwvd.  Window_length must be odd to ensure
 *         that the xwvd is not calculated between integer time
 *         samples.  If window_length is longer than or equal to the
 *         signal length, then the (normal) xwvd is generated.  If the
 *         window_length is shorter, then the psuedo-xwvd is generated
 *         (by definition).
 *
 *         NOTE: this function does not check for even window_length.
 *         It is the callers responsibility to comply.
 *
 *     unsigned window_r2
 *
 *         The smallest power of two equal to or larger than
 *         window_length.  See below for rationale and calculation.
 *
 *     unsigned window_order
 *
 *         Log2(window_r2).  See below for rationale and calculation.
 *
 * Comments:
 *
 *     Refer to comments in wvd(), above.
 *
 * Returns:
 *
 *     0       no errors
 *
 **/
int xwvd( double *signal1_real, double *signal1_imag, double
	 *signal2_real, double *signal2_imag, int signal_length, double
	 *result_real, double *result_imag, int num_slices, int
	 time_res, int window_length, int window_r2, int 
	 window_order);
