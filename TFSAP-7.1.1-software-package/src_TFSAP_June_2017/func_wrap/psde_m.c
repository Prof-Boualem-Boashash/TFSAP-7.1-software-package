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
* gateway routine for psde.c.
*************************************************/

#include <math.h>
#include "mex.h"
#include "matrix.h"
#include <string.h>

#include "arithm.h"                  /* routines for complex data */
#include "tlocal.h"                  /* local function prototypes */ 
#include "tfsa_c.h"
#include "psde.h"

void mexFunction(int nlhs, MATRIX *plhs[], int nrhs, NEED_CONST MATRIX *prhs[])
{
  double *result; 		/* output psd				*/ 
  double *signal_r;		/* real part of input signal		*/
  double *signal_i;     	/* imaginary part of input signal	*/
  
  complex *signal;
  
  int signal_length;	        /* Length of the input signal 	*/
  int window_type; 		/* Window type RECT, HANN, HAMM, BART 	*/
  int seg_len;			/* Length of the data segments 		*/
  int overlap;			/* Number of samples windows overlap	*/
  int fft_len;			/* Length of fft 			*/
  int fft_r2;                  /* Radix 2 fft length and psd length 	*/     
  int M, N, i;
  
  int n;
  char *p;     
  
  /* basic input--output number check */   
   
   if( nrhs < 5 ) {
      tfsaErr( "psde", "Not enough input arguments" );
      winNTcheck( nlhs, plhs ); 
      return;
    }
   
   if( nrhs > 5 ) {
     tfsaErr( "psde", "Too many input arguments" );
     winNTcheck( nlhs, plhs ); 
     return;
   }
   
   if( nlhs > 1 ) {
     tfsaErr( "psde", "Too many output arguments" );
     winNTcheck( nlhs, plhs ); 
     return;
   } 

   
   /* Check first input */
   
   if( !mxIsNumeric( prhs[0] ) || !MXFULL( prhs[0] ) ) {
     tfsaErr( "psde", "Input must be a vector" );
     winNTcheck( nlhs, plhs );
     return;
   }
   
   /* Get input matrix dimensions */
   
   M = ( int )mxGetM( prhs[0] );
   N = ( int )mxGetN( prhs[0] );
   
   /* Check input dimensions: A vector of unity length is considered
    * invalid; matrices are invalid. */
   
   if( (M==1 && N==1) || (M!=1 && N!=1) ) {
     tfsaErr( "psde", "Input must be a vector" );
     winNTcheck( nlhs, plhs );
     return;
   }
   
   signal_length = (int)(M>N?M:N);
   
   
   /* Check second input */
   
   if( !GoodScalar( (MATRIX *)prhs[1] ) ) {
     tfsaErr( "psde", "Segment window length must be a scalar" );
     winNTcheck( nlhs, plhs );
     return;
   }
   else      
     seg_len = (int)*(mxGetPr( prhs[1] ) );  /* Get value */
   
   
   /* Parameter value check */
   if(  seg_len < 1 ) {
     tfsaErr( "psde", "Segment length must be greater than zero" );
     winNTcheck( nlhs, plhs );
     return;
   }   
   
   /* Check window length against signal length */
   
   if( seg_len > signal_length )  {
     tfsaWarning( "psde", "Window length has been truncated to signal length" );
     seg_len = (int)signal_length;
   }

   
   /* Check third input */
   
   if( !GoodScalar( (MATRIX *)prhs[2] ) ) {
     tfsaErr( "psde", "fft length must be a scalar" );
     winNTcheck( nlhs, plhs );
     return;
   }
   else      
     fft_len = (int)*(mxGetPr( prhs[2] ) );  /* Get value */
   
   
   /* Parameter value check */
   if(  fft_len < 1 ) {
     tfsaErr( "psde", "fft length must be greater than zero" );
     winNTcheck( nlhs, plhs );
     return;
   }   
   

   fft_r2 = 2;
   /* Work out length of output signal */
   while (fft_r2 < fft_len)
     fft_r2 *= 2;

   /* check fourth input parameter */

   /* forth input */

   if( !GoodScalar( (MATRIX *)prhs[3] ) )  {
     tfsaErr( "psde", "Overlap must be a scalar" );
     winNTcheck( nlhs, plhs );
     return;
   }

   overlap = (int) *(mxGetPr(prhs[3]));
   
   /* Parameter value check */
   if(  overlap < 0 ) {
     tfsaErr( "psde", "Overlap length must be non-negative" );
     winNTcheck( nlhs, plhs );
     return;
   }   

   
   /* fifth input */
   if( MXSTRING( prhs[4] ) ) {
     
     /* Get string from Matlab following p2-51 in Mex manual */
     
     n = (int)mxGetN( prhs[4] ) + 1;  /* string length + 1 for NULL */
     p = mxCalloc( n, sizeof( char ) );
     if( p == NULL ) {
       tfsaErr( "psde", "Internal memory allocation failure" );
       winNTcheck( nlhs, plhs );
       return;
     }
    
     if( mxGetString( prhs[4], p, n ) ) {
       tfsaErr( "psde", "Could not get window type" );
       winNTcheck( nlhs, plhs );
       return;
     }
     
     
     if( !strcmp(p, "rect") || !strcmp(p, "RECT") )
       window_type = 1;         
     else if( !strcmp(p, "hann") || !strcmp(p, "HANN") )
       window_type = 2;
     else if( !strcmp(p, "hamm") || !strcmp(p, "HAMM") )
       window_type = 3;
     else if( !strcmp(p, "bart") || !strcmp(p, "BART") )
       window_type = 4;
     else {
       tfsaErr("psde", "Unknown window type");
       winNTcheck( nlhs, plhs );
       return;
     }
   }
   else 
     {
       tfsaErr("psde","Window type must be a string" );
       winNTcheck( nlhs, plhs );
       return;
     }
          
   
 /* Transferred from psde.c --- paramter checking */
   
   /* Check segment is ok */ 
   if ( (seg_len <= 1) || (seg_len > fft_len) ) {
     tfsaErr("psde", "Segment length must be > 1 and <= fft length");
     winNTcheck( nlhs, plhs );
     return;
   }
   
   /* Check overlap is ok  */
   if ( (overlap >= seg_len) || (overlap < 0) ) {
     tfsaErr("psde", "Overlap must greater than zero and less than segment length");
     winNTcheck( nlhs, plhs );
     return;
   }
   
   
   /* Dereference input */
   
   signal_r = mxGetPr(prhs[0]);	 	 /* Real part of input signal */
   if (mxIsComplex(prhs[0])) 
     signal_i = mxGetPi(prhs[0]);       /* Imag part of input signal */
   else
     signal_i = NULL;
   

   /* Create a complex signal for the input */
   signal = (complex *) mxCalloc (signal_length, sizeof(complex));
   if (signal == NULL) {
     tfsaErr("psde", "Memory allocation failed");
     winNTcheck( nlhs, plhs );
     return;
   }
   for (i=0; i < signal_length; i++)	/* signal = signal_r+j*signal_i */
     signal[i].re = signal_r[i];
   if (mxIsComplex(prhs[0])) 
     for (i=0; i<signal_length;i++)
       signal[i].im = signal_i[i];
   

      
  /* Create output vector space */
  plhs[0] = MXCREATEFULL(fft_r2, 1, sizeof(double));
  if (plhs[0] == NULL) {
    tfsaErr("psde", "Memory allocation failed");
    winNTcheck( nlhs, plhs );
    return;
  }
  result = mxGetPr(plhs[0]);  
                                                    
  psde(result, signal, signal_length, seg_len, fft_r2, overlap, window_type);
   
 }
