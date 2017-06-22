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
*
* Fast-Fourier	Transform routine
* 
*************************************************/

#include <stdlib.h>
#include <math.h>





#include "arithm.h"
#include "fft.h"

#define	TFSA_PI	3.14159265358979323846


#define	EPS 1e-7			/* arbitrary, but small	offset */

#ifdef DLLMEX							/* Add offset for PC version only... */
#define	COS(x)	cos( (double)(x	+ EPS) )
#define	SIN(x)	sin( (double)(x	+ EPS) )
#else
#define	COS(x) cos(x)
#define	SIN(x) sin(x)
#endif





/********************
 *
 * void	FFT(X, M, d)
 *
 * Calculates FFT.  Input data array must be of	length 2**M.
 *
 */
void FFT(X, M, d)
complex X[];		/* Input vector	of data	*/
int M;		/* M = log N, where N is length	of X data */
int d;		/* direction of	transform: 1=direct; -1:inverse	*/
/**************************************
 *   F F T based on FFT	in Kay's book
 **************************************/
{
    complex *Y;		 /* auxiliary array, will be created later */
    int N;		 /* number of input data */
    int num;
    int i,j,k,L;
    int LDFT, NDFT;
    complex W, SAVE;
    float ARG;
    int NP,NQ;
    double Pi = 3.14159265358979;



    N = (int) (pow(2.,(float) M)	+ 0.5);	/* calculate N=2**M */

    Y = (complex*) calloc (N, sizeof(complex));
    if (Y==NULL)
    {
        return;
    }

    /**** PUT DATA IN BIT REVERSED ORDER	******/
    Y[0]	= cmplx(X[0].re, X[0].im);
    for (i=1; i<N; i++)
    {
        num = i;
        j	= 0;
        L	= N;
        for (k = 1; k <= M; k++)
        {
            L = L/2;
            if (num >= L)
            {
                j =	j + (int) (pow (2.,(float)(k-1)) + 0.5); /* j =	j + 2**(k-1) */
                num	= num -	L;
            }
        }
        Y[i] = cmplx( X[j].re, X[j].im );
    }

    /************* Begin	computation of FFT ***********/
    LDFT	= 1;
    NDFT	= N;
    for (k = 1; k <= M; k++)
    {
        LDFT = LDFT*2;
        NDFT = NDFT/2;
        for (i = 0; i<NDFT; i++)
        {
            for (j	= 0; j < LDFT/2; j++)
            {
                /********* Determine twiddle factor	*******/
                ARG	= (float)-(2.* Pi * d *j)/LDFT;
                W =	cmplx (COS(ARG), SIN(ARG));
                /********* Determine indices for nex butterfly to be computed **/
                NP = j + LDFT*i;
                NQ = NP + LDFT/2;
                /********* Compute butterfly for the ith stage K *****/
                SAVE = add ( Y[NP],	multpl(	W, Y[NQ]) );
                Y[NQ] = sub	( Y[NP], multpl( W, Y[NQ]) );
                Y[NP] = cmplx(SAVE.re, SAVE.im);
            }
        }
    }

    for (i = 0; i < N; i++)
    {
        if (d == 1)
            X[i] =	Y[i];
        else
            X[i] =	cmult( Y[i], 1./(float)N );
    }


    free((void*)Y);
//    free(Y);

}


void compute_dft(complex* input, complex *output, int signal_length)
{
    int k, t;
    double Pi = 3.14159265358979;


    for (k = 0; k < signal_length; k++)    /* For each output element */
    {
        double sumreal = 0;
        double sumimag = 0;
        for (t = 0; t < signal_length; t++)    /* For each input element */
        {
            double angle = 2 * Pi * t * k / signal_length;
            sumreal +=  input[t].re * cos(angle) + input[t].im * sin(angle);
            sumimag += -input[t].re * sin(angle) + input[t].im * cos(angle);
        }
        output[k].re = sumreal;
        output[k].im = sumimag;
    }
}



void compute_idft(complex* input, complex *output, int signal_length)
{
    int k, t;
    double Pi = 3.14159265358979;


    for (k = 0; k < signal_length; k++)    /* For each output element */
    {
        double sumreal = 0;
        double sumimag = 0;
        for (t = 0; t < signal_length; t++)    /* For each input element */
        {
            double angle = 2 * Pi * t * k / signal_length;
            sumreal += input[t].re * cos(angle) - input[t].im * sin(angle);
            sumimag += input[t].re * sin(angle) + input[t].im * cos(angle);
        }
        output[k].re = sumreal/signal_length;
        output[k].im = sumimag/signal_length;
    }

}



void swap(complex *v1, complex *v2)
{
    complex tmp = *v1;
    *v1 = *v2;
    *v2 = tmp;
}


void fftshift(complex *data, int count)
{
    int k = 0;
    int c = (int) floor((float)count/2);
    // For odd and for even numbers of element use different algorithm
    if (count % 2 == 0)
    {
        for (k = 0; k < c; k++)
            swap(&data[k], &data[k+c]);
    }
    else
    {
        complex tmp = data[0];
        for (k = 0; k < c; k++)
        {
            data[k] = data[c + k + 1];
            data[c + k + 1] = data[k + 1];
        }
        data[c] = tmp;
    }
}

void ifftshift(complex *data, int count)
{
    int k = 0;
    int c = (int) floor((float)count/2);
    if (count % 2 == 0)
    {
        for (k = 0; k < c; k++)
            swap(&data[k], &data[k+c]);
    }
    else
    {
        complex tmp = data[count - 1];
        for (k = c-1; k >= 0; k--)
        {
            data[c + k + 1] = data[k];
            data[k] = data[c + k];
        }
        data[c] = tmp;
    }
}

