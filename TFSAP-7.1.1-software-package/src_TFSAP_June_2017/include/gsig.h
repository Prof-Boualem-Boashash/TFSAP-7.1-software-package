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
* gsig
*
* Generation of various test signals
*
*
* 
* 
*/

void linfm( double *result_real, double *result_imag, double f1, double f2, int n);

/****
 * 
 * quadfm - quadratic fm signal 
 *
 **/

void quadfm( double *result_real, double *result_imag, double f1, double f2, int n);

/****
 *
 * hypfm - hyperbolic fm signal 
 *
 **/

void hypfm( double *result_real, double *result_imag, double f1, double f2, int n);

/****
 *
 * cubicfm - cubic fm signal
 *
 **/

void cubicfm( double *result_real, double *result_imag, double f1, double f2, int n);

/****
 *
 * stepfm
 *
 * Produces a signal containing discretely stepped frequency variations.  At this time, the
 * steps are constrained to be linearly increasing, with uniform step size.  The number of 
 * steps can be specified as the last parameter.
 *
 * Note that as the number of steps increases, strange things happen.  This is because 
 * the signal begins to approximate a linear FM, but with a frequency variation of f1 to
 * f1 + 2 * (f2 - f1).
 *
 * NOTE: num_bins must be >=2, else NaN's will result.
 *
 **/
void stepfm( double *result_real, double *result_imag, double f1, double f2, int n, 
	    int num_bins);

#define WINLEN  3     /* smoothing window length for stepfm */

/****
 * sin law fm signal
 **/

 
void sinfm( double *result_real, double *result_imag, double fc, double fm, int n, double fdev );

