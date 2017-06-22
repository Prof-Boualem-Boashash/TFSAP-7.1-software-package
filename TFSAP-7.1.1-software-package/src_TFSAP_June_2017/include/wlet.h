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
 * Wavelet Transform
 *
 * Performs the forward or inverse wavelet transform using a pyramid algorithm.
 * 
 *
 */

/* output types */
#define ONE_D   1
#define TWO_D   2

/**************************************************
 *
 * int wave( double *data, double *result,  unsigned long n, int direction, unsigned num_coeff)
 *
 * Performs the Daubechies wavelet transform on data.
 *
 * double *data     Array of input data.
 *
 * double *result    Output storage.
 *
 * unsigned long n  Length of the array "data".  n must be a power of 2.
 *
 * int direction    One of:   1 ... forward transform
 *                           -1 ... inverse transform
 *
 * unsigned  num_coeff    Number of coefficients used in the definitiono of the wavelet.
 *                  Currently only values of 4, 12 and 20 are available.
 *
 * return value     is negative if an error occured, otherwise shows the number of
 *                  filtering stages required in the pyramid algorithm.
 **/
int wave( double *data, double *result, unsigned long n, int direction, unsigned num_coeff);

/* init_wave - support routine, used by wave */
int init_wave( int num_coeff, int direction);

/* wave_step - support routine, iteratively called by wave */
int wave_step( double *data, double *result, unsigned long n);

/********************
 *
 * int form_ts( double *input, double *output, unsigned long n, int m)
 *
 * Forms a time-scale matrix from the output of the function wave.
 * The effect can be seen as:
 *
 * Output from wave:
 *
 *      1 2 3 4 5 6 7 8
 *
 * After form_ts
 *
 *      1 5 7
 *      2 5 7
 *      3 6 7
 *      4 6 7
 *
 * The depth of the matrix (3 in this example) is given by the parameter m.
 *
 * Parameters:
 *
 * unsigned long n    length of input signal
 *
 * unsigned m         width of ts matrix.  This will normally equal the number of
 *                    stages in the time-scale decomposition pyramid.
 *                    The size of the ts matrix in the other direction
 *                    will be equal to n/2.  m must be less than or equal to log2(n).
 * 
 * */

int form_ts( double *input, double *output, int n); 

