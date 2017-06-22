/*************************************************
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
 * Function for computing the power spectrum density of an input signal
 * using the periodogram direct method.
 *
 * Refs:
 * ch. 5, sec. 5.7.3,
 * S. Lawrence Marple Jr.,
 * Digital Spectral Analysis with Applications.
 * Prentice-Hall Inc., Englewood Cliffs, New Jersey.
*************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mex.h"


#include "arithm.h"                  /* routines for complex data */
#include "tlocal.h"                  /* local function prototypes */ 


#include "fft.h"
#include "psde.h"
#include "tfsa_c.h"
#include "window.h"


void psde (double *result, complex *input, int signal_length, 
           int seg_len, int fft_len, int overlap, int window_type)
/*
 *  pre: fft_len must be radix 2.
 *
 */          
{
   int i, j, k;			/* Counters                           */
   complex *Y;			/* Segment minus mean and zero padded */
   complex *periodogram;	/* Pointer to the periodogram 	      */ 
   double *window;		/* pointer to the window 	      */
   double mean1=0;		/* Mean of the real part of input     */
   double mean2=0;  		/* Mean of the imag part of input     */
   int fft_order;		/* order of the fft		      */ 
   int seg_r2;			/* same as fft length                 */
   int signal_zeropad;		/* Zero padded length of signal	      */
   int num_segs;		/* Number of segments 	  	      */
   int offset;
      

   /* Calculate zero padded length of singal */
   num_segs = (int) ceil((double)signal_length/(seg_len-overlap));

   signal_zeropad = num_segs*(seg_len-overlap);
   signal_zeropad += overlap;
                                            
   /* Calculate the order of the fft for the FFT function */                                         
   fft_order = 1;       
   seg_r2 = 2;
   while (seg_r2 < fft_len) {
      fft_order++;
      seg_r2 *= 2; 
   }              
   
      
   /* Allocate space for signal minus the mean */
   Y = (complex *) mxCalloc (signal_zeropad, sizeof(complex));
   if (Y == NULL) {
     tfsaErr ("psde", "Memory allocation failed");
   }
      
       
       
   /* Copy input into zero padded array */
   for (i = 0; i < signal_length; i++) {
      Y[i].re = input[i].re;
      Y[i].im = input[i].im;
   }
  
   /* ensure that signal is zero padded */
   for (i = signal_length; i < signal_zeropad; i++) {
      Y[i].re = 0;
      Y[i].im = 0;
   }
   
   /* Allocate space for the window and the periodogram */
      window  = (double *) mxCalloc (seg_len, sizeof(double));
      periodogram  = (complex *) mxCalloc (fft_len, sizeof(complex));
   
   switch( window_type) {
      case RECT: for (k = 0; k < seg_len; k++)
 	               window[k] = 0.5; 
                 break;
      case HANN: hann(window, seg_len);
                 break;
      case HAMM: hamming(window, seg_len);
                 break;
      case BART: bartlett(window, seg_len);
                 break;
      default  : for (k=0; k < seg_len; i++)
 	               window[k] = 0.5;
   }
           
    
   periodogram[0].re = 0.3;   
                                          
   /*
    * For each data segment: apply the data window,
    *                        compute the periodogram,
    *                        sum corresponding samples.
    *  
    */

   offset = 0; 

   for (i = 0; i < num_segs; i++) {      
    
      /* Debug: Check array bounds */ 
      if ((offset+seg_len) > signal_zeropad) {
         return;
      }
        
      /* calculate mean of this segment */
      for (j = 0; j < seg_len; j++) {
         mean1 += Y[offset+j].re;
         mean2 += Y[offset+j].im;
      }
                      
      mean1 /= (double) seg_len;
      mean2 /= (double) seg_len;
      
      /* apply data window to each segment */
      for(j = 0; j < seg_len; j++) {
         periodogram[j].re = Y[offset+j].re*window[j] - mean1;
         periodogram[j].im = Y[offset+j].im*window[j] - mean2;
      } 
                                  
      /* ensure that rest of periodogram is 0 */
      for(j = seg_len; j < fft_len; j++) {
         periodogram[j].re = 0;
         periodogram[j].im = 0;
      } 
                                  
      offset += (seg_len - overlap); 

      FFT(periodogram, fft_order, 1); 

      for (j = 0; j < fft_len; j++) {
         periodogram[j].re = periodogram[j].re * periodogram[j].re +
	                         periodogram[j].im * periodogram[j].im;
       	 periodogram[j].re /= (double) fft_len;                 
   	     /* Add current segment periodogram to other segments */
	     result[j] += (double) periodogram[j].re;
      }  
   }
   /* Average periodograms */
    for (j = 0; j < fft_len; j++) {
      result[j] /= (double) num_segs;  
   }                        


}

