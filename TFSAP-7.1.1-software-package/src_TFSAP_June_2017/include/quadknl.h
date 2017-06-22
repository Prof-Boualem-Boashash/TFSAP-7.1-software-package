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
* quadknl.h
*
* Time-frequency Kernel Generation for use with Cohen's class of
* time-frequency distributions.
*
* Parameters common to all routines are described here for brevity.
*
* Parameters:
*
*     double **kernel
*
*         Points to allocated memory where kernel values will be
*         written.  kernel must be allocated as a contiguous chunk of
*         memory, and is indexed as kernel[time][lag].
*
*         Memory for kernel should be allocated by the calling function using the following code:
*
*         >  kernel = (double **)calloc( size, sizeof(double *));
*         >  kernel[0] = (double *)calloc( size * size, sizeof(double));
*         >    for (i = 1; i < size; i++)
*         >      kernel[i] = kernel[0] + size * i;
*
*         and deallocated using:
* 
*         >  free(kernel[0]);
*         >  free(kernel);
*
*     int full
*
*         Boolean value.  
*
*             1   The full kernel is generated and kernel is indexed to
*                 kernel[window_length-1][window_length-1].
*
*             0   Only the 1st quadrant is generated, and kernel is
*                 indexed to kernel[window_length/2][window_length/2].
*
*     unsigned window_length
*
*         Size of kernel.  window_length must be odd.
*
* 
**/

/* define distribution types */
#define WVD 1
#define SMOOTHED 2
#define STFT 3
#define RM 4
#define CW 5
#define BJC 6
#define ZAM 7
#define B   8
#define MB 9
#define EMB 10
#define CKD 11
#define USER 99

/* Window types */
#include "window.h"

/* prototypes */

/********************
 *
 * int wvd_kernel( double **kernel, unsigned window_length, int full)
 *
 * Generates the kernel for the (windowed) Wigner-Ville Distribution.
 *
 * Parameters:
 *
 *     double **kernel, unsigned window_length, int full
 *
 *         Described at the top of this file.
 *
 **/
int wvd_kernel( double **kernel, unsigned window_length, int full);

/********************
 *
 * int smoothedwvd_kernel( double **kernel, unsigned window_length, int
 *		       full, unsigned smooth_win_width, int smooth_win_type)
 *
 * Generates the smoothed Wigner-Ville Kernel function.
 *
 * Parameters:
 *
 *     double **kernel, unsigned window_length, int full
 *
 *         Described at the top of this file.
 *
 *     unsigned smooth_win_width
 *
 *         Width of the smoothing window.  This value must be odd.
 *
 *     int smooth_win_type
 *
 *         The type of smoothing window.  Valid types are defined at the top of this file.
 *
 **/
int smoothedwvd_kernel( double **kernel, unsigned window_length, int
		       full, unsigned smooth_win_width, int smooth_win_type);

/********************
 *
 * int stft_kernel( double **kernel, unsigned window_length, int full, int smooth_win_type)
 *
 * Generates the Spectrogram (magnitude squared of the STFT) Kernel function.
 *
 * Parameters:
 *
 *     double **kernel, unsigned window_length, int full
 *
 *         Described at the top of this file.
 *
 *     int smooth_win_type
 *
 *         The type of smoothing window.  Valid types are defined at the top of this file.
 *
 **/
int stft_kernel( double **kernel, unsigned window_length, int full,
		 unsigned smooth_win_width, int smooth_win_type);

/********************
 *
 * int rm_kernel( double **kernel, unsigned window_length, int full)
 *
 * Generates the Rihaczek-Margenau Kernel function.
 * The actual function is defined as:
 *
 *     0.5 * [ d(n+m) + d(n-m) ]
 *
 * where d() is the impulse function (d(0) = 1; d is zero otherwise).
 *
 * Parameters:
 *
 *     double **kernel, unsigned window_length, int full
 *
 *         Described at the top of this file.
 *
 **/
int rm_kernel( double **kernel, unsigned window_length, int full);

/********************
 *
 * int complex_rm_kernel( double **kernel, unsigned window_length)
 *
 * Generates the Rihaczek-Margenau Kernel function for complex distribution.
 * The actual function is defined as:
 *
 *     d(n-m)
 *
 * where d() is the impulse function (d(0) = 1; d is zero otherwise).
 *
 * Parameters:
 *
 *     double **kernel, unsigned window_length
 *
 *         Described at the top of this file.
 *
 **/
int complex_rm_kernel( double **kernel, unsigned window_length);

/********************
 *
 * int cw_kernel( double **kernel, unsigned window_length, int full, double sigma)
 *
 * Generates the Choi-Williams Kernel function.
 *
 * Parameters:
 *
 *     double **kernel, unsigned window_length, int full
 *
 *         Described at the top of this file.
 *
 *     double sigma
 *
 *         The Choi-Williams smoothing parameter.
 *
 **/
int cw_kernel( double **kernel, unsigned window_length, int full, double sigma);

/********************
 *
 * int bjc_kernel( double **kernel, unsigned window_length, int full)
 *
 * Generates the Born-Jordan-Cohen Kernel function.
 * The function used is:
 *
 *     1 / (|m| * 2 + 1)          |m| >= |n|
 *
 *     0                          otherwise
 *
 * Parameters:
 *
 *     double **kernel, unsigned window_length, int full
 *
 *         Described at the top of this file.
 **/
int bjc_kernel( double **kernel, unsigned window_length, int full);

/********************
 *
 * int zam_kernel( double **kernel, unsigned window_length, int full, int a)
 *
 * Generates the Zhao-Atlas-Marks Kernel function.
 * The function used is:
 *
 *
 * Parameters:
 *
 *     double **kernel, unsigned window_length, int full
 *
 *         Described at the top of this file.
 **/
int zam_kernel( double **kernel, unsigned window_length, int full, double a);

/********************
 *
 * int b_kernel( double **kernel, unsigned window_length, int full, double beta)
 *
 * Generates the B Kernel function.
 *
 * Parameters:
 *
 *     double **kernel, unsigned window_length, int full
 *
 *         Described at the top of this file.
 *
 *     double beta
 *
 *         The B kernel smoothing parameter.
 *
 **/
int b_kernel( double **kernel, unsigned window_length, int full, double beta);

/********************
 *
 * int mb_kernel( double **kernel, unsigned window_length, int full, double alpha)
 *
 * Generates the modified B Kernel function.
 *
 * Parameters:
 *
 *     double **kernel, unsigned window_length, int full
 *
 *         Described at the top of this file.
 *
 *     double alpha
 *
 *         The modified B kernel smoothing parameter.
 *
 **/
int mb_kernel( double **kernel, unsigned window_length, int full, double alpha);



/********************
 *
 * int emb_kernel( double **kernel, unsigned window_length, int full, double alpha)
 *
 * Generates the Extended modified B Kernel function.
 *
 * Parameters:
 *
 *     double **kernel, unsigned window_length, int full
 *
 *         Described at the top of this file.
 *
 *     double alpha & beta
 *
 *         The Extended modified B kernel smoothing parameter.
 *
 **/
int emb_kernel( double **kernel, unsigned window_length, int full, double alpha, double beta);


























