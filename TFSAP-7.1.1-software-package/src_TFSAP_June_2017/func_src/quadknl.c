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
* The following references 2 should be cited whenever this script is used:
* [1] B. Boashash, Samir Ouelha, Designing time-frequency  and time-scale
* features for efficient classification of non stationary signals: a tutorial
* review with performance comparison, Digital Signal Processing, In Press.
* [2] B. Boashash, Samir Ouelha, Efficient Software Matlab Package for the calculation
* of Time-Frequency Distributions related Time-Scale methods and the extraction of
* signal characteristics, SoftwareX, In Press.
* In addition, the following 3rd reference is relevant for use of the executable
* [3] B. Boashash (ed.), Time-Frequency Signal Analysis and Processing, 2nd Edition
* (London: Elsevier / Academic Press, December 2015); ISBN 978-0-12-398499-9. Book website:
* http://booksite.elsevier.com/9780123984999/
*
* Description:
*
* Time-frequency Kernel Generation for	use with quadratic
* time-frequency distributions.
*************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "mex.h"

#include "fft.h"
#include "arithm.h"		     /*	for complex vectors */
#include "tlocal.h"		     /*	local function prototypes */

#include "quadknl.h"
#include "tfsa_c.h"



/********************/
int wvd_kernel(	double **kernel, unsigned window_length, int full)
{
    unsigned m;


    if (full)
    {
        memset( kernel[0], 0, (size_t)(sizeof(double)*(unsigned long)window_length*(unsigned long)window_length));
        for	(m = 0;	m < window_length; m++)
            kernel[0][m] = 1;
    }
    else
    {
        memset( kernel[0], 0, (size_t)(sizeof(double)*((unsigned long)window_length+1)*((unsigned long)window_length+1)/4));
        for	(m = 0;	m < (window_length+1)/2; m++)
            kernel[0][m] = 1;
    }

    return 0;
}

/********************/
int smoothedwvd_kernel(	double **kernel, unsigned window_length, int
                        full, unsigned smooth_win_width,	int smooth_win_type)
{
    unsigned m;
    unsigned n;
    double *window;
    int i;
    int wc;			/* window centre */



    window  = (double *)mxCalloc(smooth_win_width, sizeof(double));

    if( window ==NULL )
    {
        tfsaErr("quadknl", "Internal memory allocation failure: window" );
        return 1;
    }

    switch( smooth_win_type )
    {
    case RECT:
        for( i=0; i<(int)smooth_win_width; i++ )
        {
            window[i]	= 1.0 /	(double)smooth_win_width;
        }
        break;

    case HANN:
        hann(window, smooth_win_width);
        break;

    case HAMM:
        hamming(window, smooth_win_width);
        break;

    case BART:
        bartlett(window, smooth_win_width);
        break;

    default:

        /* This should never be executed...	parameters checked in gateway. */

        for	(i=1; i<(int)smooth_win_width; i++)
        {
            window[i]	= 1.0 /	(double)smooth_win_width;
        }
        break;
    }


    wc = (int)smooth_win_width/2;

    if (full)
    {
        memset( kernel[0], 0, (size_t)(sizeof(double)*(unsigned long)window_length*(unsigned long)window_length));

        for	( m = 0; m < window_length; m++)
        {
            kernel[0][m] = window[wc];
            for ( n =	1; n < (smooth_win_width + 1)/2; n++)
            {
                kernel[n][m] = window[wc+n];
                kernel[window_length-n][m] = window[wc-n];
            }
        }

    }
    else
    {
        memset( kernel[0], 0, (size_t)(sizeof(double)*((unsigned long)window_length+1)*((unsigned long)window_length+1)/4));

        for	( m = 0; m < (window_length+1)/2; m++)
            for ( n =	0; n < (smooth_win_width + 1)/2; n++)
                kernel[n][m] = window[wc+n];
    }


    mxFree( (void	*)window);

    return 0;
}



int stft_kernel( double	**kernel,
                 unsigned window_length,
                 int full,
                 unsigned smooth_win_width,
                 int smooth_win_type)
{
    unsigned m;
    unsigned n;
    unsigned max;
    double *window;
    int wc;
    int i;


    max =	(smooth_win_width-1)/2;


    window  = (double*)mxCalloc(window_length, sizeof(double));

    if( window ==	NULL )
    {
        tfsaErr("quadknl", "Memory allocation	failed:	window");
        return 1;
    }

    if (full)
    {

        memset( kernel[0], 0, (size_t)(sizeof(double)*(unsigned long)window_length*(unsigned long)window_length));


        /* smooth_win_width	must be	odd */
        /* if s_w_w	= 11, then centre is index [5] */
        wc = (int)smooth_win_width/2;

        switch( smooth_win_type )
        {
        case RECT:
            for (i=0; i<(int)smooth_win_width; i++)
                window[i] = 1.0	/ (double)smooth_win_width;
            break;
        case HANN:
            hann(window, smooth_win_width);
            break;
        case HAMM:
            hann(window, smooth_win_width);
            /* TMP !!!! */
            /* hamming(window, smooth_win_width); */
            break;
        case BART:
            bartlett(window, smooth_win_width);
            break;
        default:

            for (i=1; i<(int)smooth_win_width; i++)
            {
                window[i] = 1.0	/ (double)smooth_win_width;
            }
            break;
        }



        kernel[0][0] = window[wc]*window[wc];
        for( n = 1;	n <= max; n++)
        {
            kernel[n][0] = window[wc+n]*window[wc+n];
            kernel[window_length - n][0] = window[wc-n]*window[wc-n];
        }


        for	(m = 1; m <= max; m++)
        {
            kernel[0][m] = window[wc+m]*window[wc-m];
            kernel[0][window_length - m] = window[wc-m]*window[wc+m];

            for(n = 1; n <= (max - m); n++)
            {
                kernel[n][m] = window[wc+n+m]*window[wc+n-m];;
                kernel[n][window_length-m] = window[wc+n-m]*window[wc+n+m];
                kernel[window_length-n][m] = window[wc-n+m]*window[wc-n-m];
                kernel[window_length-n][window_length-m] = window[wc-n-m]*window[wc-n+m];;
            }
        }


    }    
    else
    {
        memset( kernel[0], 0, (size_t)(sizeof(double)*((unsigned long)window_length+1)*((unsigned long)window_length+1)/4));


        /* smooth_win_width must be odd	*/
        /* if s_w_w = 11, then centre is index [5] */
        wc = (int)smooth_win_width/2;

        switch(smooth_win_type)
        {
        case RECT:
            for (i=0; i<(int)smooth_win_width; i++)
                window[i] = 1.0 / (double)smooth_win_width;
            break;
        case HANN:
            hann(window, smooth_win_width);
            break;
        case HAMM:
            hamming(window, smooth_win_width);
            break;
        case BART:
            bartlett(window, smooth_win_width);
            break;
        default:

            for (i=1; i<(int)smooth_win_width; i++)
            {
                window[i] = 1.0 / (double)smooth_win_width;
            }
            break;
        }

        /* Changed from kernel[n][m] = window[wc+n]; */
        /* 20/10/03. JOT */

        for ( m = 0; m <= max; m++)
        {
            for ( n = 0; n <= (max-m); n++)
                kernel[n][m] = window[wc+n+m]*window[wc+n-m];
        }

    }

    mxFree( (void	*)window);

    return 0;
}

/********************/
int rm_kernel( double **kernel,	unsigned window_length,	int full)
{
    unsigned i;

    if (full)
    {
        memset( kernel[0], 0, (size_t)(sizeof(double)*(unsigned long)window_length*(unsigned long)window_length));
        kernel[0][0] = 1;
        for	(i = 1;	i<(window_length+1)/2; i++)
        {
            kernel[i][i] = 0.5;
            kernel[window_length-i][i] = 0.5;
            kernel[i][window_length-i] = 0.5;
            kernel[window_length-i][window_length-i] = 0.5;
        }
    }
    else
    {
        memset( kernel[0], 0, (size_t)(sizeof(double)*((unsigned long)window_length+1)*((unsigned long)window_length+1)/4));
        kernel[0][0] = 1;
        for	(i = 1;	i<(window_length+1)/2; i++)
            kernel[i][i] = 0.5;
    }

    return 0;
}

/********************/
int complex_rm_kernel( double **kernel,	unsigned window_length)
{
    unsigned i;

    memset( kernel[0], 0,	(size_t)(sizeof(double)*(unsigned long)window_length*(unsigned long)window_length));
    kernel[0][0] = 1;
    for (i = 1; i<(window_length+1)/2; i++)
    {
        kernel[i][i] = 1;
        kernel[window_length-i][window_length-i] = 1;
    }
    return 0;
}

/********************/
int cw_kernel( double **kernel,	unsigned window_length,	int full, double sigma)
{
    int m, n;
    double x = sqrt( sigma / 3.14159265358979323846) / 2;
    double y = -sigma / 4;
    int size = (window_length+1)/2;

    if (full)
    {

        kernel[0][0] = 1;

        /* m = 0 */
        for( n = 1;	n < size; n++)
        {
            kernel[n][0] = 0;
            kernel[window_length - n][0] = 0;
        }
        /* n = 0 */
        for( m = 1;	m < size; m++)
        {
            kernel[0][m] = x / m;
            kernel[0][window_length -	m] = x / m;
        }

        for	(m = 1;	m < size; m++)
            for ( n =	1; n < size; n++)
            {
                kernel[n][m] = x / m * exp( y *	n * n /	m / m);
                kernel[window_length-n][m] = kernel[n][m];
                kernel[n][window_length-m] = kernel[n][m];
                kernel[window_length-n][window_length-m] = kernel[n][m];
            }
    }
    else
    {
        kernel[0][0] = 1;
        for (m = 1; m < size; m++)
            for ( n	= 0; n < size; n++)
                kernel[n][m] = x / m * exp( y	* n * n	/ m / m);
    }

    return 0;
}

/********************/
int bjc_kernel(	double **kernel, unsigned window_length, int full)
{
    int m, n;
    int size = (window_length+1)/2;

    if (full)
    {
        memset( kernel[0], 0, (size_t)(sizeof(double)*(unsigned long)window_length*(unsigned long)window_length));
        kernel[0][0] = 1;

        /* kernel =	0 for all m = 0, n != 0	*/

        /* n = 0 */
        for( m = 1;	m < size; m++)
        {
            kernel[0][m] = 1.0 / (double)((m*2)+1);
            kernel[0][window_length -	m] = kernel[0][m];
        }

        for	(m = 1;	m < size; m++)
            for ( n =	1; n <=	m; n++)
            {
                kernel[n][m] = 1.0 / (double)((m*2)+1);
                kernel[window_length-n][m] = kernel[n][m];
                kernel[n][window_length-m] = kernel[n][m];
                kernel[window_length-n][window_length-m] = kernel[n][m];
            }
    }
    else
    {
        memset( kernel[0], 0, (size_t)(sizeof(double)*((unsigned long)window_length+1)*((unsigned long)window_length+1)/4));
        kernel[0][0] = 1;
        for (m = 1; m < size; m++)
            for ( n	= 0; n <= m; n++)
                kernel[n][m] = 1.0 / (double)((m*2)+1);
    }

    return 0;
}

/********************/
int zam_kernel(	double **kernel, unsigned window_length, int full, double a)
{
    int m, n;
    int size = (window_length+1)/2;
    double alpha = log(100.0)/(2.0 * size	* size);
    double wt;



    if (full)
    {
        memset( kernel[0], 0, (size_t)(sizeof(double)*(unsigned long)window_length*(unsigned long)window_length));
        kernel[0][0] = 1;

        /* kernel =	0 for all m = 0, n != 0	*/

        /* n = 0 */
        for( m = 1;	m < size; m++)
        {
            kernel[0][m] = exp( -2.0 * alpha * m * m);
            kernel[0][window_length -	m] = kernel[0][m];
        }

        for	(m = 1;	m < size; m++)
        {
            wt = exp(-2.0 * alpha * m	* m);
            for ( n =	1; n <=	(int)(m/a); n++)
            {
                kernel[n][m] = wt;
                kernel[window_length-n][m] = kernel[n][m];
                kernel[n][window_length-m] = kernel[n][m];
                kernel[window_length-n][window_length-m] = kernel[n][m];
            }
        }
    }
    else
    {
        memset( kernel[0], 0, (size_t)(sizeof(double)*((unsigned long)window_length+1)*((unsigned long)window_length+1)/4));
        kernel[0][0] = 1;
        for (m = 1; m < size; m++)
        {
            wt = exp(-2.0 *	alpha *	m * m);
            for ( n	= 0; n <= (int)(m/a); n++)
                kernel[n][m] = wt;
        }
    }

    return 0;
}

/********************/
int b_kernel( double **kernel, unsigned window_length, int full, double beta)
{

    int  m, n;
    int size = (window_length + 1) / 2;

    if (full)
    {
        memset( kernel[0], 0, (size_t)(sizeof(double)*(unsigned long)window_length*(unsigned long)window_length));
        kernel[0][0] = 0;

        // m = 0 /
        for (n = 1; n < size; n++)
        {
            kernel[n][0] = 0;
            kernel[window_length - n][0] = 0;
        }

        // n = 0 /
        for (m = 1; m < size; m++)
        {
            kernel[0][m] = pow(m, beta);
            kernel[0][window_length - m] = pow(m, beta);
        }

        for	(m = 1;	m < size; m++)
            for (n = 1; n < size; n++)
            {
                kernel[n][m] = pow((m/pow(cosh(n), 2)), beta);
                kernel[window_length - n][m] = kernel[n][m];
                kernel[n][window_length - m] = kernel[n][m];
                kernel[window_length - n][window_length - m] = kernel[n][m];
            }
    }
    else
    {
        memset( kernel[0], 0, (size_t)(sizeof(double)*((unsigned long)window_length+1)*((unsigned long)window_length+1)/4));
        for (m = 0; m < size; m++)
            for (n = 0; n < size; n++)
                kernel[n][m] = pow((m/pow(cosh(n), 2)), beta);
    }

    return 0;
}


/********************/
int mb_kernel( double **kernel, unsigned window_length, int full, double alpha)
{
    int m, n;
    int size = (window_length + 1) / 2;
    double tmp_gammaln, k_alpha;


    /* Function lgamma is not standard <math.h> ANSI C. */
    /* For the momemnt will use MATLABs one till I get a */
    /* proper free 'free' function.                */


    MATRIX *gamma_lhs, *gamma_rhs;

    gamma_rhs = MXCREATEFULL (1, 1, REAL);
    *( mxGetPr( gamma_rhs ) ) = 2 * alpha;

    if( mexCallMATLAB(1, &gamma_lhs, 1, &gamma_rhs, "gammaln" ) )
    {
        tfsaErr("quadknl", "Unable to call gammaln MATLAB function" );
        return 1;
    }

    tmp_gammaln = (double)*( mxGetPr( gamma_lhs ) );
    mxDestroyArray( (MATRIX *)gamma_rhs );
    mxDestroyArray( (MATRIX *)gamma_lhs );


    k_alpha = exp( tmp_gammaln )/(pow(2,(2*alpha-1))*pow(exp( tmp_gammaln ),2));


    if (full)
    {
        memset( kernel[0], 0, (size_t)(sizeof(double)*(unsigned long)window_length*(unsigned long)window_length));
        kernel[0][0] = k_alpha;

        /* m = 0 */
        for (n = 1; n < size; n++)
        {
            kernel[n][0] = k_alpha/pow(cosh(n),(2*alpha));
            kernel[window_length - n][0] = kernel[n][0];
        }

        /* n = 0 */
        for (m = 1; m < size; m++)
        {
            kernel[0][m] = k_alpha;
            kernel[0][window_length - m] = kernel[0][m];
        }

        for	(m = 1;	m < size; m++)
            for (n = 1; n < size; n++)
            {
                kernel[n][m] = k_alpha/pow(cosh(n),(2*alpha));
                kernel[window_length - n][m] = kernel[n][m];
                kernel[n][window_length - m] = kernel[n][m];
                kernel[window_length - n][window_length - m] = kernel[n][m];
            }
    }
    else
    {


        memset( kernel[0], 0, (size_t)(sizeof(double)*((unsigned long)window_length+1)*((unsigned long)window_length+1)/4));
        for (m = 0; m < size; m++)
            for (n = 0; n < size; n++)
            {
                kernel[n][m] = k_alpha/pow(cosh(n),(2*alpha));
            }
    }

    return 0;
}

int emb_kernel( double **kernel, unsigned window_length, int full, double alpha, double beta)
{
    int  i,j,m, n;
    int size = (window_length + 1) / 2;

    double *kernel_time;
    double *kernel_lag;

    double *data1, *data2;
    mxArray *gamma_lhs, *gamma_rhs;

    // Create an array for the output's data /
    data2 = (double *) mxMalloc(2 * sizeof(double));

    data2[0]=beta;
    data2[1]=window_length;

    gamma_rhs = mxCreateDoubleMatrix (1, 2, mxREAL);
    mxSetPr(gamma_rhs, data2);

    if( mexCallMATLAB(1, &gamma_lhs, 1, &gamma_rhs, "gamma_wraper" ) )
    {
        tfsaErr("quadknl", "Unable to call gamma_wrapper MATLAB function");
        return 1;
    }

    // Get the data passed in /
    data1 = mxGetPr(gamma_lhs);

    kernel_time = (double *)calloc(size, sizeof(double));
    if (kernel_time==NULL)
    {
        printf("\n Unable to allocate Memory");
        return 0;
    }


    kernel_lag = (double *) mxCalloc (size, sizeof(double));

    if (kernel_lag==NULL)
    {
        printf("\n Unable to allocate Memory");
        return 0;
    }
    // cosh function calculation for time
    for (m = 0; m < size; m++)
    {
        kernel_time[m] = 1/pow(cosh(m),(2*alpha));
    }

    for	(m = 0; m < size; m++)
    {
        kernel_lag[m] = data1[m];
    }

    // EMBD Kernel


    if (full)
    {
        memset( kernel[0], 0, (size_t)(sizeof(double)*(unsigned long)window_length*(unsigned long)window_length));

        for (m = 0; m < size; m++)
        {
            for (n = 0; n < size; n++)
            {
                kernel[n][m] = kernel_time[n]*kernel_lag[m];
                if (m>0)
                {
                    kernel[n][window_length-m] = kernel[n][m];

                    if (n>0)
                    {
                        kernel[window_length-n][m-1] = kernel[n][m];
                        kernel[window_length-n][window_length-m] = kernel[n][m];
                    }
                }
            }
        }
    }
    else
    {
        memset( kernel[0], 0, (size_t)(sizeof(double)*((unsigned long)window_length+1)*((unsigned long)window_length+1)/4));

        for (m = 0; m < size; m++)
        {
            for (n = 0; n < size; n++)
            {
                kernel[n][m] = kernel_time[n]*kernel_lag[m];
            }
        }
    }

    return 0;
}




//***********************************














