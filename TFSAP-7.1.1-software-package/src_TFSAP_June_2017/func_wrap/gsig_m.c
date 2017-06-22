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
*
* Description:
*	
* gateway routine for gsig.c.
*************************************************/

#include <string.h>
#include <math.h>
#include "mex.h"
#include "arithm.h"		     /*	routines for complex data */
#include "tlocal.h"		     /* local function prototypes */
#include "tfsa_c.h"                  /* needed for preprocessor macros */
#include "gsig.h"


void mexFunction(int nlhs, MATRIX *plhs[], int nrhs, NEED_CONST MATRIX *prhs[])
{
  double f1;				/* normalise frequency lower bound */
  double f2;				/* normalise frequency upper bound */

  double cf;				/* central frequency */
  double mf;				/* modulation frequency	*/
  double fdev;				/* frequency deviation */

  int num;
  short	sig_type;			/* 1 = real, otherwise complex */
  double *result_r, *result_i;

  int num_bins;

  int n;
  char *p;

  /* basic input--output check */

  if( nrhs < 5 ) {
    tfsaErr( "gsig", "Not enough input arguments" );
    winNTcheck(	nlhs, plhs );
    return;
  }

  if( nrhs > 6 ) {     /* Six or more inputs */
    tfsaErr( "gsig", "Too many input arguments"	);
    winNTcheck(	nlhs, plhs );
    return;
  }

  if( nlhs > 1 ) {
    tfsaErr( "gsig", "Too many output arguments" );
    winNTcheck( nlhs, plhs );
    return;
  }

  /* First input: Check if string */

  if( !MXSTRING(prhs[0]) ) {
    tfsaErr( "gsig", "Data type must be a string" );
    winNTcheck( nlhs, plhs );
    return;
  }

  /* Get string from Matlab env. See p2-51 in Mex manual */

  n = (int)mxGetN( prhs[0] ) + 1;  /* string length + 1 for NULL */
  p = mxCalloc( n, sizeof( char ) );
  if( p == NULL ) {
    tfsaErr( "gsig", "Internal memory allocation failure" );
    winNTcheck( nlhs, plhs );
    return;
  }

  if( mxGetString( prhs[0], p, n ) ) {
    tfsaErr( "gsig", "Could not get data type" );
    winNTcheck( nlhs, plhs );
    return;
  }

  /*  Check to see if valid data type given.  Respond accordingly */

  /***************************** Data type 1 *********************/

  if( !strcmp( p, "lin" )  || !strcmp( p, "LIN" )   ||
     !strcmp( p, "quad")   || !strcmp( p, "QUAD" )  ||
     !strcmp( p, "cubic" ) || !strcmp( p, "CUBIC" ) ||
     !strcmp( p, "hyp" )   || !strcmp( p, "HYP" ) ) {


  if( nrhs > 5 ) {     /* Must have exactly five inputs */
    tfsaErr( "gsig", "Too many input arguments" );
    winNTcheck( nlhs, plhs );
    return;
  }



    /* Check second input parameter: f1 */

    if( !GoodScalar( (MATRIX *)prhs[1] ) ) {
      tfsaErr( "gsig", "Second input parameter must be a scalar" );
      winNTcheck( nlhs, plhs );
      return;
    }
    else 		   /* Get value -- real part */
      f1  = mxGetScalar( prhs[1] );

    /* Check third input parameter: f2 */

    if( !GoodScalar( (MATRIX *)prhs[2] ) ) {
      tfsaErr( "gsig", "Second input parameter must be a scalar" );
      winNTcheck( nlhs, plhs );
      return;
    }
    else			   /* Get value -- real part */
      f2  = mxGetScalar( prhs[2] );

    /* Check fourth input parameter: num_samples */

    if( !GoodScalar( (MATRIX *)prhs[3] ) ) {
      tfsaErr( "gsig", "Fourth input parameter must be a scalar" );
      winNTcheck( nlhs, plhs );
      return;
    }
    else
      num = (int) mxGetScalar( prhs[3] );

    /* Check value: num */
    if( num < 1 ) {
      tfsaErr( "gsig", "Number of samples must be greater than zero" );
      winNTcheck( nlhs, plhs );
      return;
    }

    /* Check fifth input parameter: sig_type */

    if( !GoodScalar( (MATRIX *)prhs[4] ) ) {
      tfsaErr( "gsig", "sig_type flag must be 0 (analytic) or 1 (real)" );
      winNTcheck( nlhs, plhs );
      return;
    }
    else
      sig_type = (short) mxGetScalar( prhs[4] );

    /* Allocate space for the output */

    if( sig_type == 1 )
      plhs[0] = MXCREATEFULL( num, 1, REAL ); /* REAL/COMPLEX  in <matrix.h> */
    else
      plhs[0] = MXCREATEFULL( num, 1, COMPLEX );

    if( plhs[0] == NULL ) {
      tfsaErr( "gsig", "Memory allocation failure" );
      winNTcheck( nlhs, plhs );
      return;
    }

    result_r = mxGetPr( plhs[0] );

    if( result_r == NULL ) {
      tfsaErr( "gsig", "result_r is NULL!" );
      winNTcheck( nlhs, plhs );
      return;
    }


    if( sig_type == 1 )			/* Real data */
      result_i = NULL;
    else
    {

      result_i = mxGetPi( plhs[0] );	/* Complex data */

    if( result_i == NULL ) {
      tfsaErr( "gsig", "result_i is NULL!" );
      winNTcheck( nlhs, plhs );
      return;
    }

}


    /* call respective data_type1 function */

    if( !strcmp (p, "lin") || !strcmp (p, "LIN" ) )
    { linfm( result_r, result_i, f1, f2, num ); }

    else
    if( !strcmp (p, "quad") || !strcmp (p, "QUAD") )
     {quadfm( result_r, result_i, f1, f2, num );}

    else if( !strcmp (p, "hyp") || !strcmp (p, "HYP") )
      {hypfm( result_r, result_i, f1, f2, num);}

    else cubicfm( result_r, result_i, f1, f2, num);

  }

  /*****************  Data_type2 ******************/

  else
    if( !strcmp( p, "sin" ) || !strcmp( p, "SIN" ) ) {

      /* Check for correct number of input arguments for sin */
      if( nrhs < 6 )
 {
	  tfsaErr( "gsig", "Not enough input arguments for type 'sin'" );
	  winNTcheck( nlhs, plhs );
	  return;
	}

      /* Check second input parameter: cf */

      if( !GoodScalar( (MATRIX *)prhs[1] ) ) {
 	tfsaErr( "gsig", "Second input parameter must be a scalar" );
	winNTcheck( nlhs, plhs );
	return;
      }
      else			     /*	Get value -- real part */
	cf  = mxGetScalar( prhs[1] );

      /* check cf value	*/


      /* Check third input parameter: mf */

      if( !GoodScalar( (MATRIX *)prhs[2] ) ) {
	tfsaErr( "gsig", "Second input parameter must be a scalar" );
	winNTcheck( nlhs, plhs );
	return;
      }
      else			     /* Get value -- real part */
	mf  = mxGetScalar( prhs[2] );

      /* Check fourth input parameter: num_samples */

      if( !GoodScalar( (MATRIX *)prhs[3] ) ) {
	tfsaErr( "gsig", "Fourth input parameter must be a scalar" );
	winNTcheck( nlhs, plhs );
	return;
      }
      else
	num = (int) mxGetScalar( prhs[3] );

      /* Check value: num */
      if( num <	1 ) {
	tfsaErr( "gsig", "Number of samples must be greater than zero" );
	winNTcheck( nlhs, plhs );
	return;
      }

      /* Check fifth input parameter: sig_type */

      if( !GoodScalar( (MATRIX *)prhs[4] ) ) {
	tfsaErr( "gsig", "sig_type flag must be 0 (analytic) or 1 (real)" );
	winNTcheck( nlhs, plhs );
	return;
      }
      else
	sig_type = (short) mxGetScalar( prhs[4] );


      /* Check sixth input parameter: fdev */

      if( !GoodScalar( (MATRIX *)prhs[5] ) ) {
	tfsaErr( "gsig", "Second input parameter must be a scalar" );
	winNTcheck( nlhs, plhs );
	return;
      }
      else			     /*	Get value -- real part */
	fdev  =	mxGetScalar( prhs[5] );

      /* Allocate space for the output */

      if( sig_type == 1 )
	plhs[0] = MXCREATEFULL(num, 1, REAL);
      else
	plhs[0] = MXCREATEFULL(num, 1, COMPLEX);

      if( plhs[0] == NULL ) {
	tfsaErr( "gsig", "Memory allocation failure" );
	winNTcheck( nlhs, plhs );
	return;
      }

      result_r = mxGetPr( plhs[0] );
      if( sig_type == 1	)			/* Real	data */
	result_i = NULL;
      else
	result_i = mxGetPi( plhs[0] );	  /* Complex data */

      /* invoke function */

      sinfm( result_r, result_i, cf, mf, num, fdev );

    }
    else if( !strcmp( p, "step" ) || !strcmp( p, "STEP" ) ) {

      /*************** data_type3 ***********************/

      /* input parameter check */

      if( nrhs < 6 ) {
	tfsaErr( "gsig", "Not enough input arguments" );
	winNTcheck( nlhs, plhs );
	return;
      }


      /* Check second input parameter: f1 */

      if( !GoodScalar( (MATRIX *)prhs[1] ) ) {
	tfsaErr( "gsig", "Second input parameter must be a scalar" );
	winNTcheck( nlhs, plhs );
	return;
      }
      else			     /*	Get value -- real part */
	f1  = mxGetScalar( prhs[1] );

      /* Check third input parameter: f2 */

      if( !GoodScalar( (MATRIX *)prhs[2] ) ) {
	tfsaErr( "gsig", "Second input parameter must be a scalar" );
	winNTcheck( nlhs, plhs );
	return;
      }
      else			     /*	Get value -- real part */
	f2  = mxGetScalar( prhs[2] );

      /* Check fourth input parameter: num_samples */

      if( !GoodScalar( (MATRIX *)prhs[3] ) ) {
	tfsaErr( "gsig", "Fourth input parameter must be a scalar" );
	winNTcheck( nlhs, plhs );
	return;
      }
      else
	num = (int) mxGetScalar( prhs[3] );

      /* Check value: num */
      if( num <	1 ) {
	tfsaErr( "gsig", "Number of samples must be greater than zero" );
	winNTcheck( nlhs, plhs );
	return;
      }

      /* Check fifth input parameter: sig_type */

      if( !GoodScalar( (MATRIX *)prhs[4] ) ) {
	tfsaErr( "gsig", "sig_type flag must be 0 (analytic) or 1 (real)" );
	winNTcheck( nlhs, plhs );
	return;
      }
      else
	sig_type = (short) mxGetScalar( prhs[4] );

      /* Check fifth input parameter: sig_type */

      if( !GoodScalar( (MATRIX *)prhs[5] ) ) {
	tfsaErr( "gsig", "Fifth input parameter must be a scalar" );
	winNTcheck( nlhs, plhs );
	return;
      }
      else
	num_bins = (int) mxGetScalar( prhs[5] );

      if( num_bins < 2 ) {
	tfsaErr( "gsig", "Number of frequency steps must be greater than or equal to two" );
	winNTcheck( nlhs, plhs );
	return;
      }


      /* Allocate space for the output */

      if( sig_type == 1	)
	plhs[0] = MXCREATEFULL(num, 1, REAL);
      else
	plhs[0] = MXCREATEFULL(num, 1, COMPLEX);

      if( plhs[0] == NULL ) {
	tfsaErr( "gsig", "Memory allocation failure" );
	winNTcheck( nlhs, plhs );
	return;
      }

      result_r = mxGetPr( plhs[0] );
      if( sig_type == 1 )			/* Real data */
	result_i = NULL;
      else
	result_i = mxGetPi( plhs[0] );   /* Complex data */

      /* Invoke function */

      stepfm (result_r, result_i, f1, f2, num, num_bins);

    }
    else {
      tfsaErr( "gsig", "Unknown data type" );
      winNTcheck( nlhs, plhs );
      return;
    }
}
