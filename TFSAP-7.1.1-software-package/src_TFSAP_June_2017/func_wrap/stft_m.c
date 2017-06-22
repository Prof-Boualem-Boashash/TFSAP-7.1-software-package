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
* gateway routine for stft.c.
*************************************************/

#include "mex.h"
#include "matrix.h"

#include "stft.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *tfd_re, *tfd_im;
	double *signal;
	double *time_vector, *frequency_vector;
    char *p;
	int signal_length, window_length, Fe, window_overlap, Nfft, stft_or_spec, M, N, nbFrames, n, window_type, nbRows;




	/******************* basic input--output number check ************************/

    if( nrhs != 7 )
        mexErrMsgTxt( "Seven inputs required." );

    if (nlhs >3)
		mexErrMsgTxt("Too many output arguments");



  /********************* Check first input *****************************************/

  /* Get input matrix dimensions */

  M = ( int )mxGetM( prhs[0] );
  N = ( int )mxGetN( prhs[0] );

  /* Check input dimensions: A vector of unity length is considered
   * invalid; matrices are invalid. */

    if( (M==1 && N==1) || (M!=1 && N!=1) )
    {
        mexErrMsgTxt( "Input must be a vector" );
        return;
    }

 signal_length = (int)(M>N?M:N);
 signal = mxGetPr( prhs[0]);

  if (mxIsComplex(prhs[0]))
  {
    mexErrMsgTxt( "Input must be a real vector" );
  }


  /***************************************** Check second input ***********************/

    M = mxGetM( prhs[1] );
    N = mxGetN( prhs[1] );


    if( (M!=1 || N!=1) )
    {
      mexErrMsgTxt( "Sampling frequency must be a scalar" );
      return;
    }
    else
    {
    Fe = (int)mxGetScalar(prhs[1]);
    }
  /************************* Check third input *************************************/

    M = mxGetM( prhs[2] );
    N = mxGetN( prhs[2] );


    if( (M!=1 || N!=1) )
    {
      mexErrMsgTxt( "Sampling frequency must be a scalar" );
      return;
    }
    else
    window_length = (int)*(mxGetPr( prhs[2] ) );  /* Get value */


  /* Parameter value check */
  if( window_length < 1 ) {
    mexErrMsgTxt("Window length must be greater than zero" );
    return;
  }

  /* Check window length against signal length */

  if( window_length > signal_length )
  {
    mexWarnMsgTxt("Window length has been truncated to signal length" );
    window_length = (int)signal_length;
  }


  /******************** Check fourth input ********************************************/

    if( mxIsChar( prhs[3] ) )
    {

      n = (int)mxGetN( prhs[3] ) + 1;  /* string length + 1 for NULL */
      p = mxCalloc( n, sizeof( char ) );
      if( p == NULL )
	{
	  mexErrMsgTxt("Internal memory allocation failure" );
	  return;
	}

      if( mxGetString( prhs[3], p, n ) )
	{
	  mexErrMsgTxt( "Could not get smoothing window type" );
	  return;
	}

      if( !strcmp( p, "rect" ) || !strcmp( p, "RECT" ) )
	window_type = RECT;
      else if( !strcmp( p, "hann" ) || !strcmp( p, "HANN" ) )
	window_type = HANN;
      else if( !strcmp( p, "hamm" ) || !strcmp( p, "HAMM" ) )
	window_type = HAMM;
      else if( !strcmp( p, "bart" ) || !strcmp( p, "BART" ) )
	window_type = BART;
      else
	{
	  mexErrMsgTxt("Unknown smoothing window type" );
	  return;
	}
    }
  else
    {
      mexErrMsgTxt("Smoothing window type must be a string" );
      return;
    }

/******************** Check fifth input ********************************************/

    M = ( int )mxGetM( prhs[4] );
    N = ( int )mxGetN( prhs[4] );

    if( (M!=1 || N!=1) )
    {
        mexErrMsgTxt( "Overlapping must be a scalar" );
        return;
    }
    else
    {
        window_overlap = (int)mxGetScalar( prhs[4]);
    }

      if( window_overlap > window_length )
    {
        mexWarnMsgTxt("Window overlap has been truncated to window length minus one" );
        window_overlap = (int)window_length-1;
    }

/******************** Check sixth input ********************************************/

    M = mxGetM( prhs[5] );
    N = mxGetN( prhs[5] );

    if( (M!=1 || N!=1) )
    {
        mexErrMsgTxt( "Nfft flag must be a scalar");
        return;
    }
    else
    {
       Nfft = (int)mxGetScalar(prhs[5]) ;
    }

    if( window_length > Nfft )
    {
        mexWarnMsgTxt("fft length has been fixed to window length minus one" );
        Nfft = (int)window_length;
    }

 /******************** Check seventh input ********************************************/

    M = mxGetM( prhs[6] );
    N = mxGetN( prhs[6] );

    if( (M!=1 || N!=1) )
    {
        mexErrMsgTxt( "STFT/Spectrogram flag must be a scalar");
        return;
    }
    else
    {
       stft_or_spec = (int)mxGetScalar(prhs[6]) ;
    }

/******************** Check first output ********************************************/

    nbRows = (Nfft/2)+1;
    nbFrames = (int)ceil((signal_length-window_overlap)/(window_length-window_overlap));


    plhs[0] = mxCreateDoubleMatrix( nbRows, nbFrames, mxCOMPLEX );
    tfd_re = mxGetPr( plhs[0] );
    tfd_im = mxGetPi( plhs[0] );



    if( !plhs[0] )
    {
        mexErrMsgTxt("Memory allocation failed" );
        return;
    }

	/******************** Check second output ********************************************/
	plhs[1] = mxCreateDoubleMatrix(nbFrames,1, mxREAL );
	time_vector = mxGetPr( plhs[1] );

	/******************** Check third output ********************************************/
	plhs[2] = mxCreateDoubleMatrix( 1, nbRows, mxREAL );
	frequency_vector = mxGetPr( plhs[2] );


	/********************* Call C-subroutine for the Short Time Fourier Transform *************/

stft(signal, tfd_re, tfd_im,time_vector, frequency_vector, signal_length, Fe, window_type, window_length, window_overlap, Nfft, stft_or_spec);
}
