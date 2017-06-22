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
* S-method
*
* This function computes the S-method time-frequency representation.
*
*************************************************/

#include "specSM.h"



void specSM(double *signal, double *tfd_real, int L, int signal_length, int window_type, int window_length, int window_overlap,int fft_length)
{
    int l,i, n, m, offset=0;
    complex temp;
    int nbRows = (fft_length/2)+1;


    float temp1 = (signal_length-window_overlap)/(window_length-window_overlap);

    int nbFrames = ceil(temp1);
    double *stftReal = calloc(nbRows*nbFrames, sizeof(double));
    double *stftImag = calloc(nbRows*nbFrames, sizeof(double));
    double *tfd_imag = calloc(nbRows*nbFrames, sizeof(double));
    double *time_vector = calloc(nbFrames, sizeof(double));
    double *frequency_vector = calloc(nbRows, sizeof(double));
    double *fenetreP = calloc(2*L+1, sizeof(double));
    complex *stftComplex = calloc(nbRows*nbFrames, sizeof(complex));



    switch(window_type)
    {
    case RECT:
        for (l=0; l<2*L+1; l++)
        {
            fenetreP[l] = 1.0/(double)(2*L+1);
        }
        break;
    case HANN:
        hann(fenetreP, 2*L+1);
        break;
    case HAMM:
        hamming(fenetreP, 2*L+1);
        break;
    case BART:
        bartlett(fenetreP, 2*L+1);
        break;
    default:
        for (l=0; l<2*L+1; l++)
        {
            fenetreP[l] = 1.0/(double)(2*L+1);
        }
        break;
    }

    stft(signal, stftReal, stftImag , time_vector, frequency_vector, signal_length, 1, window_type, window_length, window_overlap, fft_length, 1);

    for (i=0; i<nbRows*nbFrames; i++)
    {
        stftComplex[i].re = stftReal[i];
        stftComplex[i].im = stftImag[i];
    }


    for (n=0; n<nbFrames; n++)
    {
        for (m=0; m<nbRows; m++)
        {
            for (l=-L; l<L+1; l++)
            {
                if ((m+l<nbRows&&m+l>=0)&&(m-l<nbRows&&m-l>=0))
                {
                    offset = m+n*nbRows;
                    temp = multpl ((stftComplex[(m+l)+n*nbRows]), (conj(stftComplex[(m-l)+n*nbRows])));
                    tfd_real[offset] = tfd_real[offset]+fenetreP[L+l]*temp.re;
                    tfd_imag[offset] = tfd_imag[offset]+fenetreP[L+l]*temp.im;
                }
            }
        }
    }


//    if (fichier != NULL)
//    {
//        for (i = 0; i<nbRows*nbFrames; i++)
//        {
//
//            fprintf (fichier, "%5lf ", tfd_real[i]);
//        }
//        fclose(fichier);
//    }

    free(stftReal);
    free(stftImag);
    free(tfd_imag);
    free(fenetreP);
    free(stftComplex);
    free(time_vector);
    free(frequency_vector);

}





