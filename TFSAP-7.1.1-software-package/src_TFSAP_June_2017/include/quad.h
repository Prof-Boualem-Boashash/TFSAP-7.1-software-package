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
* quadtfd
*
* Quadratic Class of TF Distributions
*
**/

/***
 *
 * int quad_complexG_asymmG( double *signal_real, double *signal_imag,
 *          unsigned signal_length, double *result_real, double
 *          *result_imag, unsigned num_slices, unsigned time_res,
 *          unsigned window_length, unsigned window_r2, unsigned
 *          window_order, double **G_real, double **G_imag)
 *
 * int quad_realG_asymmG( double *signal_real, double *signal_imag,
 *          unsigned signal_length, double *result_real, unsigned
 *          num_slices, unsigned time_res, unsigned window_length,
 *          unsigned window_r2, unsigned window_order, double **G_real)
 *
 * int quad_complexG_symmG( double *signal_real, double *signal_imag,
 *          unsigned signal_length, double *result_real, double
 *          *result_imag, unsigned num_slices, unsigned time_res,
 *          unsigned window_length, unsigned window_r2, unsigned
 *          window_order, double **G_real, double **G_imag)
 *
 * int quad_realG_symmG( double *signal_real, double *signal_imag,
 *          unsigned signal_length, double *result_real, unsigned
 *          num_slices, unsigned time_res, unsigned window_length,
 *          unsigned window_r2, unsigned window_order, double **G_real)
 * 
 * Description:
 * 
 *     Computes the tf distribution of the input signal. An analytic
 *     signal generator is called if the input signal is real.  Four
 *     versions of the program appear, so that faster algorithms can
 *     be used when certain conditions on the kernel G are met.
 *
 *     When G is real, a version with _realG in the function name
 *     should be used.  The _realG functions can assume the resulting
 *     distribution will be real and need perform only half the number
 *     of ffts otherwise required.
 *
 *     If G is symmetric in both time and lag, use the _symmG
 *     variants.  For these functions, only the first quadrant of G
 *     need be specified in the kernel array passed.
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
 *     unsigned signal_length
 *
 *         Length of signal_real and _imag arrays.
 *
 *     double *result_real
 *     double *result_imag
 *
 *         Should point to allocated memory to store real and
 *         imaginary components of the computed distribution.  If the
 *         distribution can be assumed to be real, then result_imag
 *         can be passed as NULL, and only the real component will be
 *         returned.  These result arrays are returned holding
 *         window_r2 * num_slices doubles each; num_slices of
 *         window_r2 consecutive frequency values (beginning at 0 sec
 *         and 0 Hz samples).
 *
 *     unsigned num_slices
 *
 *         The number of time instants at which the distribution will
 *         be computed.  This is closely related to time_res.  See
 *         below for more details and how to calculate it.
 *
 *     unsigned time_res
 *
 *         The time-spacing between successive frequency slices of the
 *         distribution.  This is closely related to num_slices.  See
 *         below for more details.
 *
 *     unsigned window_length
 *
 *         window_length is the size of the determining function G,
 *         where G is defined from -window_length/2 to window_length/2
 *         in both dimensions.  window_length must be odd.
 *
 *         NOTE: this function does not check for odd window_length.
 *         It is the callers responsibility to comply.
 *
 *     unsigned window_r2
 *
 *         The smallest power of two equal to or larger than
 *         window_length.  See below for rationale and calculation.
 *         Extra zero-padding at the fft stage of the analysis may be
 *         specified using a larger than normal value here.
 * 
 *     unsigned window_order
 *
 *         Log2(window_r2).  See below for rationale and calculation.
 *
 *     complex **G_real, **G_imag
 *
 *         The kernel function defining the particular member of
 *         Cohen's class.  G is defined in the time-lag plane, and
 *         indexing occurs as G[n][m] where n is the time variable and
 *         m the lag variable.  Negative values of n and m are mapped
 *         as G[window_length-n][window_length-m].  If G is symmetric,
 *         the variation with _symmG in the function title may be
 *         called.  These functions will not index beyond
 *         G[(window_length+1)/2][(window_length+1)/2], as they
 *         retrieve first quadrant information only and deduce the
 *         other quadrants by symmetry.
 *         
 *         See below in Comments for code to allocate memory for G.
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
 *     following code.  Refer to comments in function for more detail.
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
 *     >       result_r = (double *)malloc( sizeof(double) * window_r2 * num_slices);
 *     >       result_i = (double *)malloc( sizeof(double) * window_r2 * num_slices);
 *     >     // end //
 *
 *     Note that result_i does not need to be allocated if one of the
 *     _realG versions of the program is being called.
 *
 *     Memory for the determining function can be allocated using the
 *     following code:
 *
 *     >     // addressing is G[n][m] //
 *     >
 *     >     if ( "using _symmG" )
 *     >       size = (window_length+1)/2;
 *     >     else
 *     >       size = window_length;
 *     >
 *     >     #ifdef MATLAB_MEX_FILE
 *     >     G_real = (double **)mxCalloc( size, sizeof(double *));
 *     >     G_real[0] = (double *)mxCalloc( size * size, sizeof(double));
 *     >     for (i = 1; i < size; i++)
 *     >       G_real[i] = G_real[0] + size * i;
 *     >     #else
 *     >     G_real = (double **)calloc( size, sizeof(double *));
 *     >     G_real[0] = (double *)calloc( size * size, sizeof(double));
 *     >     for (i = 1; i < size; i++)
 *     >       G_real[i] = G_real[0] + size * i;
 *     >     #endif
 *
 *     As an example, the wvd can be generated using a determining
 *     function defined as:
 *
 *     >     memset( G_real[0], 0, sizeof(double)*size*size);
 *     >     memset( G_imag[0], 0, sizeof(double)*size*size);
 *     >     for (i = 0; i < size; i++)
 *     >       G_real[0][i] = 1;
 *
 *     Deallocate G using:
 *
 *     >     #ifdef MATLAB_MEX_FILE
 *     >       mxFree( G[0]);
 *     >       mxFree( G);
 *     >     #else
 *     >       free(G[0]);
 *     >       free(G);
 *     >     #endif
 *
 * Returns:
 *
 *     0       no errors
 *     1       memory allocation error
 *
 **/

int quad_complexG_asymmG( 
		  double *signal_real, double *signal_imag,
		  unsigned signal_length, 
		  double *result_real, double *result_imag,
		  unsigned num_slices, unsigned time_res, 
		  unsigned window_length, unsigned window_r2, unsigned window_order,
		  double **G_real, double **G_imag);

int quad_realG_asymmG(
		  double *signal_real, double *signal_imag,
		  unsigned signal_length, 
		  double *result_real,
		  unsigned num_slices, unsigned time_res, 
		  unsigned window_length, unsigned window_r2, unsigned window_order,
		  double **G_real);

int quad_complexG_symmG(
		  double *signal_real, double *signal_imag,
		  unsigned signal_length, 
		  double *result_real, double *result_imag,
		  unsigned num_slices, unsigned time_res, 
		  unsigned window_length, unsigned window_r2, unsigned window_order,
		  double **G_real, double **G_imag);

int quad_realG_symmG(
		  double *signal_real, double *signal_imag,
		  unsigned signal_length, 
		  double *result_real,
		  unsigned num_slices, unsigned time_res, 
		  unsigned window_length, unsigned window_r2, unsigned window_order,
		  double **G_real);









