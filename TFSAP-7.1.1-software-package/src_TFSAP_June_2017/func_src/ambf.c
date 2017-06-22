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
* [2] B. Boashash, Samir Ouelha, Efficient Software Matlab Package for the calculation of Time-Frequency Distributions related
* Time-Scale methods and the extraction of signal characteristics, SoftwareX, 2017.
* In addition, the following 3rd reference is relevant for use of the executable
* [3] B. Boashash (ed.), Time-Frequency Signal Analysis and Processing, 2nd Edition
* (London: Elsevier / Academic Press, December 2015); ISBN 978-0-12-398499-9. Book website:
* http://booksite.elsevier.com/9780123984999/
*
* Description:
*
* This is the latest version of the C code using
* MEX-files to compute the ambiguity function
* of a signal.
*
* This file takes into account the lag-resolution
* and the windowing so that a faster computation is
* obtained. See Amb.m for more details.
*
*************************************************/
#include <math.h>
#include <ctype.h>
#include <stdlib.h>


#include "mex.h"


#include "arithm.h"		     /*	routines for complex data */
#include "tlocal.h"		     /*	local function prototypes */


//#include "ambf.h"
#include "analyt.h"
#include "fft.h"
#include "tfsa_c.h"

#define pi 3.1415926535

/********************/
int Amb(double *signal_real, double *signal_imag, int signal_length,
        int fft_length, int r1, double *result_real, double *result_imag, int lag_res,
        int wind_length)
{
    int i;
    complex *sig;
    complex *Kernel1, *Kernel2, *Kernel11, *Kernel22;
    int m, n, index1, index2;

    /* read in input signal into complex array */


    sig =	(complex *) mxCalloc (signal_length, sizeof(complex));
    if (sig==NULL)
    {
        tfsaErr("Amb", "Memory allocation failed");
        return 1;
    }


    if (signal_imag == NULL)
    {
        for	(i = 0;	i < signal_length; i++)
        {
            sig[i].re	= signal_real[i];
            sig[i].im	= 0.0;
        }
        default_sigana(sig,signal_length);
    }
    else
        for	(i=0; i	< signal_length; i++)
        {
            sig[i].re	= signal_real[i];
            sig[i].im	= signal_imag[i];
        }

    /* Here begins the real computation */

    /* Memory allocation for the kernel vectors  */

    Kernel1=(complex*) mxCalloc(fft_length,sizeof(complex));
    Kernel11=(complex*) mxCalloc(fft_length,sizeof(complex));
    Kernel2=(complex*) mxCalloc(fft_length,sizeof(complex));
    Kernel22=(complex*) mxCalloc(fft_length,sizeof(complex));

    for(m=0; m<(wind_length-1)/2; m=lag_res+m) /* m is used for the lag */
    {

        for(n=0; n<signal_length; n++)
        {
            index1=n+m;
            index2=n-m;
            if (index1<signal_length)
            {
                if (index2>=0)
                {
                    Kernel1[n].re=sig[index1].re*sig[index2].re+sig[index1].im*sig[index2].im;
                    Kernel1[n].im=-sig[index1].re*sig[index2].im+sig[index1].im*sig[index2].re;
                    Kernel2[n].re=Kernel1[n].re;
                    Kernel2[n].im=-Kernel1[n].im;
                }
            }
        }


        /* FFT of the kernel for each lag m, one FFT is for the positive lag
         * and the other FFT is for the corresponding negative lag (-m),
         * we use K(n,-m)=K*(n,m) */

        FFT(Kernel1,r1,1);
        FFT(Kernel2,r1,1);

        /* fftshift of the vectors Kernel1 an Kernel2 */
        for(i=0; i<fft_length/2; i++)
        {
            Kernel11[fft_length/2+i].re=Kernel1[i].re;
            Kernel11[fft_length/2+i].im=Kernel1[i].im;
            Kernel22[i].re=Kernel2[fft_length/2+i].re;
            Kernel22[i].im=Kernel2[fft_length/2+i].im;
        }
        for(i=0; i<fft_length/2; i++)
        {
            Kernel11[i].re=Kernel1[fft_length/2+i].re;
            Kernel11[i].im=Kernel1[fft_length/2+i].im;
            Kernel22[fft_length/2+i].re=Kernel2[i].re;
            Kernel22[fft_length/2+i].im=Kernel2[i].im;
        }

        /* Building the AF by rows and re-initialising the
        * the kernel vectors for the next FFT use. */

        for(i=0; i<fft_length; i++)
        {
            *(result_real+(-m+(wind_length+1)/2)+i*wind_length)=Kernel22[i].re;
            *(result_imag+(-m+(wind_length+1)/2)+i*wind_length)=Kernel22[i].im;
            *(result_real+(m+(wind_length+1)/2)+i*wind_length)=Kernel11[i].re;
            *(result_imag+(m+(wind_length+1)/2)+i*wind_length)=Kernel11[i].im;
        }

        for(i=0; i<fft_length; i++)
        {
            Kernel1[i].re=0;
            Kernel1[i].im=0;
            Kernel2[i].re=0;
            Kernel2[i].im=0;
        }
    }


    return 0;
}





