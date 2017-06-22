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
*
* Description:
*
* Gateway Routine for ambf.c
*************************************************/
#include <math.h>
#include "mex.h"
#include "matrix.h"

#include "arithm.h"		     /*	routines for complex data */
#include "tlocal.h"		     /*	local function prototypes */

#include "ambf.h"
#include "tfsa_c.h"


void mexFunction( int nlhs, MATRIX *plhs[], int nrhs, NEED_CONST MATRIX *prhs[] )
{
    double *result_r, *result_i; 	      /* pointer to results matrix data */
    double *signal_r, *signal_i;
    int signal_length;
    int M,N,r1,r2;
    int fft_length;
    int lag_res, wind_length;

    /* basic input--output number	check */

    if( nrhs != 3)
    {
        tfsaErr( "Amb", "Three input arguments must exist" );
        winNTcheck( nlhs, plhs );
        return;
    }

    if( nlhs >1 )
    {
        tfsaErr( "Amb", "Only one output argument must exist." );
        winNTcheck(	nlhs, plhs );
        return;
    }

    /* Check first input */

    if( !mxIsNumeric( prhs[0] ) || !MXFULL( prhs[0] ) )
    {
        tfsaErr( "Amb", "Input signal must	be a vector" );
        winNTcheck(	nlhs, plhs );
        return;
    }

    /* Get input matrix dimensions */

    M = (	int )mxGetM( prhs[0] );
    N = (	int )mxGetN( prhs[0] );

    /* Check input dimensions: A vector of unity length is considered
     * invalid; matrices are invalid. */

    if( (M==1 && N==1) ||	(M!=1 && N!=1) )
    {
        tfsaErr( "Amb", "Input signal must	be a vector" );
        winNTcheck(	nlhs, plhs );
        return;
    }

    signal_length	= (unsigned)(M>N?M:N);

    /* Get lag resolution value */


    if( !GoodScalar( (MATRIX *)prhs[1] ))
    {
        tfsaErr( "Amb", "Lag resolution must be an integer");
        winNTcheck( nlhs, plhs );
        return;
    }
    else
        lag_res= (int) *(mxGetPr( prhs[1] ) );

    /* lag resolution check */
    if( lag_res<1)
    {
        tfsaErr("Amb","lag resolution must be greater than zero");
        winNTcheck(nlhs, plhs);
        return;
    }

    /* Get window length value */

    if( !GoodScalar( (MATRIX *)prhs[2] ))
    {
        tfsaErr( "Amb", "window length must be an integer");
        winNTcheck( nlhs, plhs );
        return;
    }
    else
        wind_length= (int) *(mxGetPr( prhs[2] ) );

    /* window length check */
    if( wind_length<1)
    {
        tfsaErr("Amb","window length  must be greater than zero");
        winNTcheck(nlhs, plhs);
        return;
    }

    /* calculate	radix-2	value equal or above (2*wind_length) */
    r1= 0;
    r2 = 1;
//   while( r2 < (2*wind_length)) {
    while( r2 < (signal_length))
    {

        r1++;
        r2<<= 1;
    }
    fft_length=r2;


    /* Allocate memory for the output variable result
     *
     * NB: The ambiguity function is complex.
     *
     */


    plhs[0] = MXCREATEFULL( wind_length, fft_length, 1 );

    result_r = mxGetPr( plhs[0] );
    result_i = mxGetPi( plhs[0] );

    /* Dereference input	*/

    signal_r = mxGetPr( prhs[0] );
    if( mxIsComplex(prhs[0]) )
        signal_i =	mxGetPi( prhs[0] );
    else
        signal_i =	NULL;

    Amb( signal_r, signal_i, signal_length, fft_length, r1, result_r, result_i, lag_res, wind_length);

}

