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
*
* Short-time Fourier transform with overlap using direct transformation following the usual definition.
*
* 
*
*************************************************/
#include "stft.h"



void stft(double *signal, double *tfd_real, double *tfd_imag , double *time_vector, double *frequency_vector, int signal_length, int Fe, int window_type, int window_length, int window_overlap,int Nfft, int stft_or_spec)
{
    int l, i, offset=0, offset2 = 0, offset3 = 0, window_order, window_r2;
    complex *temp = calloc(Nfft , sizeof(complex));
    double *window = calloc(window_length, sizeof(double));
    double sommeFenetre = 0;
    int nbRows = (Nfft/2)+1;
    int nbFrames = (int)ceil((signal_length-window_overlap)/(window_length-window_overlap));

    double *table_temp = calloc(Nfft*nbFrames, sizeof(double));

    for (l=0;l< nbFrames;l++)
    {
        time_vector[l] = (double)l*((double)window_length-(double)window_overlap)/(double)Fe+((double)window_length-1)/(2*(double)Fe);
    }

    for (l=0;l< nbRows;l++)
    {
        frequency_vector[l] = (double)l*(double)Fe/(double)Nfft;
    }



     switch(window_type)
     {
         case RECT:
             for (l=0;l<window_length;l++)
             {
                 window[l] = 1.0/(double)window_length;
             }
             break;
         case HANN:
             hann(window, window_length);
             break;
         case HAMM:
             hamming(window, window_length);
             break;
         case BART:
             bartlett(window, window_length);
             break;
         default:
            for (l=0;l<window_length;l++)
            {
                 window[l] = 1.0/(double)window_length;
            }
             break;
     }

     sommeFenetre = sommeTableau(window, window_length);

    for(i=0;i<window_length;i++)
     {
        for(l=0;l<nbFrames;l++)
        {
            offset = i+l*Nfft;
            table_temp[offset] = signal[(window_length-window_overlap)*l+i]*window[i];
        }
     }



    window_order = 0;
    window_r2 = 1;
    while( window_r2 < Nfft)
    {
        window_order++;
        window_r2 <<= 1;
    }


    for(l=0;l<nbFrames;l++)
    {
        for (i=0;i<Nfft;i++)
        {
            offset2 = i+l*Nfft;
            temp[i].re= table_temp[offset2];
            temp[i].im=0;
        }

       FFT(temp, window_order, 1);

        for (i=0;i< nbRows;i++)
        {
            offset3 = i+l*nbRows;
            tfd_real[offset3]= temp[i].re;
            tfd_imag[offset3]= temp[i].im;
        }
    }

        if (!stft_or_spec)
        {
            for (i=0;i<nbFrames*nbRows;i++)
            {
                tfd_real[i] = tfd_real[i] * tfd_real[i] + tfd_imag[i] * tfd_imag[i];
                tfd_real[i] = tfd_real[i]/sommeFenetre;
                tfd_imag[i] = 0;
            }

        }

    free(temp);
    free(window);
    free(table_temp);

}



double sommeTableau(double tableau[], int tailleTableau)
{
    int i;
    double somme=0;

    for (i=0; i<tailleTableau; i++)
    {
        somme = somme + tableau[i];
    }

    return somme;
}



