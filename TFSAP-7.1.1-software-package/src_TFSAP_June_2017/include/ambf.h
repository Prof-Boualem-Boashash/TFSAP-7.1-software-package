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
 * Amb
 *
 * Program for generating ambiguity function.
 *
 * If MATLAB_MEX_FILE is defined, then the Matlab memory allocation
 * routines are called instead of the system calls.  This is usefull
 * so that the memory is released if the user presses ^C during
 * calculation.
 * 
 **/

#ifndef ARITHM_H
#include "arithm.h"
#endif

/****
 *
 * int Amb( double *signal_real, double *signal_imag, int
 *          signal_length, int fft_length, int r1, double *result_real, double *result_imag)
 *
 * Description:
 * 
 *     Computes the ambiguity function of the input signal.  An analytic signal
 *     generator is called if the input signal is real. 
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
 *     double *result_real
 *     double *result_imag
 *
 *         Should point to allocated memory to store the computed real part
 *         and imaginary part of the ambiguity function.
 *         These result arrays are returned holding signal_length*fft_length
 *         doubles each.
 *         
 ***/



int Amb(double *signal_real, double *signal_imag, int signal_length,
	int fft_length, int r1, double *result_real, double *result_imag, int lag_res, int wind_length);

