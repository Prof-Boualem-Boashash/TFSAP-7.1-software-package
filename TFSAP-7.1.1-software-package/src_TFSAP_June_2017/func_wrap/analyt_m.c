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
* analyt.c -- gateway routine for analytic signal generation
*************************************************/
#include <stdio.h>
#include "mex.h"

#include "arithm.h"		     /*	routines for complex data */
#include "tlocal.h"		     /*	local function prototypes */
#include "analyt.h"		     /*	prototype for default_sigana() */
#include "tfsa_c.h"                  /* needed for preprocessor defs */

void mexFunction( int nlhs, MATRIX *plhs[], int	nrhs, NEED_CONST MATRIX *prhs[] )
{

  double *signal_r;
  complex *signal;
  unsigned signal_length;

  double *result_r, *result_i;

  int M, N;

  /* Check input arguments */
  switch( nrhs ) {

  case(	0 ):   /* No input */

    tfsaErr( "analyt", "Not enough input arguments" );
    winNTcheck(	nlhs, plhs ); return;
    break;

  case(	1 ):   /* One input */

    /* Check if	input is well defined */

    if(	!mxIsNumeric( prhs[0] )	|| mxIsComplex(	 prhs[0] )
       || !MXFULL( prhs[0] ) ) {
      tfsaErr( "analyt", "Input must be a real vector" );
      winNTcheck( nlhs,	plhs );
      return;
    }

    /* Get input matrix	dimensions */

    M =	( int )mxGetM( prhs[0] );
    N =	( int )mxGetN( prhs[0] );

    /* Check input dimensions: A vector	of unity length	is considered
       invalid;	matrices are invalid. */


    if(	(M==1 && N==1) || (M!=1	&& N!=1) ) {
      tfsaErr( "analyt", "Input must be a vector" );
      winNTcheck( nlhs,	plhs );
      return;
    }

    /* Choose the larger row or	column dimension */

    signal_length = (unsigned)(	M>N ? M:N );


    /* Dereference input: Real part only */
    signal_r = mxGetPr(	prhs[0]	);

    /* printf( "signal_r[0]=%f\n", signal_r[0] ); */
    break;

  default:   /*	More than one input */

    tfsaErr( "analyt", "Too many input arguments" );
    winNTcheck(	nlhs, plhs );
    return;
    break;
  }


  /* Check output arguments.  All input	parameters are correct at this
     stage. */

  switch( nlhs ) {

  case(	0 ):
  case(	1 ):

    /* allocate	output row matrix, and check memory allocation */

    plhs[0]=MXCREATEFULL( signal_length, 1, COMPLEX );

    if(	plhs[0]	== NULL	) {
      tfsaErr( "analyt", "Internal memory allocation failure" );
      winNTcheck( nlhs,	plhs );
      return;
    }

    /* Set <result> to plhs[0].	NB: If nlhs == 0, then result is
     * returned	in <ans> (Matlab takes care of this) after exiting the
     * gateway.	*/

    result_r = mxGetPr(	plhs[0]	);
    result_i = mxGetPi(	plhs[0]	);

    break;

  default:  /* more than one output */

    tfsaErr( "analyt", "Too many output	arguments" );
    winNTcheck(	nlhs, plhs );
    return;
    break;
  }

  /* Convert input to type complex */
  signal = (complex *)mxCalloc(	signal_length, sizeof( complex ) );

  if( signal ==	NULL )	{
    tfsaErr( "analyt", "Internal memory	allocation failure" );
    winNTcheck(	nlhs, plhs );
    return;
  }

  RealToComplex( signal, signal_r, signal_length );

  
  default_sigana( signal, signal_length	);

  FromComplex( signal, result_r, result_i, signal_length );

}


