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
*
* In addition, the following 3rd reference is relevant for use of the executable
* [3] B. Boashash (ed.), Time-Frequency Signal Analysis and Processing, 2nd Edition
* (London: Elsevier / Academic Press, December 2015); ISBN 978-0-12-398499-9. Book website:
* http://booksite.elsevier.com/9780123984999/
*
* Description:
*
*Function for Generating Time-Frequency Distribution based
on Compact Support Estimation
*************************************************/

#include <string.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "arithm.h"		     /* routines for complex data */
#include "tlocal.h"		     /* local function prototypes */
#include "tfsa_c.h"
#include "cmpt.h"



// Function for Generating Time-Frequency Distribution based
// on Compact Support Estimation
void cmpt_filtering( double *signal_real,
        int signal_length, int window_r2,
        double **G_real, double *result_r) {

    int i, j, f1, f2, t1, a, b, ind, ind_i, ind_j;
    int N, size_window;
    double *z_real, *z_imag;

	// declare variable R for one slice of TFD //
    complex *R;

    double *result_real;
    double *result_imag;

    double *result_real_invt;
    double *result_imag_invt;

    double *res_real_ifft_invt;
    double *res_imag_ifft_invt;

    double *res_final_real;
    double *res_final_imag;

    double **K_real;
    double **K_imag;

    double *K_real_fft_time_1_D;
    double *K_imag_fft_time_1_D;

    double *K_real_inv_fft_time_1_D;
    double *K_imag_inv_fft_time_1_D;

    double *real_final_time_1_D;

    MATRIX *fft_lhs, *fft_rhs;
    MATRIX *fft_lhs_2, *fft_rhs_2;
    MATRIX *fft_lhs_final, *fft_rhs_final;

    double real_z_a, imag_z_a;
    double real_z_b_con, imag_z_b_con;

    /* allocate storage for one slice of TFD */
    R =	(complex *) mxCalloc(window_r2, sizeof(complex));
    if (R==NULL) {
        mexPrintf( "\nMemory allocation	failed\n"	);
    }

    result_real = (double *) mxMalloc( window_r2 * signal_length * sizeof(double));
    result_imag = (double *) mxMalloc( window_r2 * signal_length * sizeof(double));

    result_real_invt = (double *) mxMalloc( signal_length * window_r2 * sizeof(double));
    result_imag_invt = (double *) mxMalloc( signal_length * window_r2 * sizeof(double));

    res_real_ifft_invt = (double *) mxMalloc( signal_length * window_r2 * sizeof(double));
    res_imag_ifft_invt = (double *) mxMalloc( signal_length * window_r2 * sizeof(double));

    res_final_real = (double *) mxMalloc( window_r2 * signal_length * sizeof(double));
    res_final_imag = (double *) mxMalloc( window_r2 * signal_length * sizeof(double));

    z_real = (double *) mxMalloc(signal_length*2 * sizeof(double));
    z_imag = (double *) mxMalloc(signal_length*2 * sizeof(double));

    fft_rhs = MXCREATEFULL( signal_length, window_r2, COMPLEX );
    fft_rhs_2 = MXCREATEFULL( signal_length, window_r2, COMPLEX );
    fft_rhs_final = MXCREATEFULL( window_r2, signal_length, COMPLEX );

    K_real =mxMalloc( window_r2*sizeof(double*));
    if(K_real==NULL) {
        printf("\nError allocating memory\n");
    }
    //*  allocate each row  */
    for(i = 0; i < window_r2; i++) {
        K_real[i] = mxMalloc( signal_length*sizeof(double));
    }
    if(K_real[i-1]==NULL) {
        printf("\nError allocating memory\n");
    }
    K_imag = mxMalloc( window_r2*sizeof(double*));
    if(K_imag==NULL) {
        printf("\nError allocating memory\n");
    }
    //*  allocate each row  */
    for(i = 0; i < window_r2; i++) {
        K_imag[i] = mxMalloc( signal_length*sizeof(double));
    }
    if(K_imag[i-1]==NULL) {
        printf("\nError allocating memory\n");
    }

    // Calling Analytic Signal Function//
    analytic_signal( signal_real, z_real, z_imag, signal_length);

    N = 2*signal_length;
    size_window = (signal_length+1)/2;


    // Calculating Instantaneous Auto-correlation Function (IAF in Time-Lag Domain) //
    for (i = 0; i <= signal_length-1; (i=i+1)) {
        t1=i;

        // Calculating IAF for each time-slice
        //R[j] =	z[t+j].z*[t-j],  with corrected indices:
        //
        for (j = 0; j <= size_window-1; j++) {
            f1=t1+j;
            f2=t1-j;

            // Computing Modulus for indexes
            a = Modulus_Index(f1, N);
            b = Modulus_Index(f2, N);

            real_z_a = z_real[a];
            imag_z_a = z_imag[a];

            real_z_b_con = z_real[b];
            imag_z_b_con = -z_imag[b];

            R[j].re = real_z_a*real_z_b_con - imag_z_a*imag_z_b_con;
            R[j].im = real_z_a*imag_z_b_con + imag_z_a*real_z_b_con;

            if (j>0) {
                R[window_r2-j].re = R[j].re;
                R[window_r2-j].im = -R[j].im;
            }
        }

        for (ind = 0; ind <= window_r2-1; ind++) {
            result_real[ind+window_r2*i]=R[ind].re;
            result_imag[ind+window_r2*i]=R[ind].im;
        }
    }

    for (ind_i=0;ind_i<window_r2;ind_i++) {
        for (ind_j=0;ind_j<signal_length;ind_j++) {

            result_real_invt[signal_length*ind_i+ind_j] = result_real[ind_i+window_r2*ind_j];
            result_imag_invt[signal_length*ind_i+ind_j] = result_imag[ind_i+window_r2*ind_j];

        }
    }

	// Take fourier transform in time domain to form Doppler-Lag Ambugity Domain //
    mxSetPr(fft_rhs, result_real_invt);
    mxSetPi(fft_rhs, result_imag_invt);
    // Calling Matlab function fft for taking fourier transform in time domain
    mexCallMATLAB(1, &fft_lhs, 1, &fft_rhs, "fft" );
    K_real_fft_time_1_D = mxGetPr(fft_lhs);
    K_imag_fft_time_1_D = mxGetPi(fft_lhs);

    // Filtering Doppler-Lag Ambugity Domain with Compact Support Kernel (G_real) //
    for (ind_i=0;ind_i<window_r2;ind_i++) {
        for (ind_j=0;ind_j<signal_length;ind_j++) {
            res_real_ifft_invt[signal_length*ind_i+ind_j] = K_real_fft_time_1_D[signal_length*ind_i+ind_j] * G_real[ind_i][ind_j];
            res_imag_ifft_invt[signal_length*ind_i+ind_j] = K_imag_fft_time_1_D[signal_length*ind_i+ind_j] * G_real[ind_i][ind_j];
        }
    }

    // Take inverse fourier transform in time domain to form time-lag Domain
	mxSetPr(fft_rhs_2, res_real_ifft_invt);
    mxSetPi(fft_rhs_2, res_imag_ifft_invt);
    // Calling Matlab function ifft for taking inverse fourier transform in time domain
    mexCallMATLAB(1, &fft_lhs_2, 1, &fft_rhs_2, "ifft" );
    K_real_inv_fft_time_1_D = mxGetPr(fft_lhs_2);
    K_imag_inv_fft_time_1_D = mxGetPi(fft_lhs_2);

    for (ind_i=0;ind_i<window_r2;ind_i++) {
        for (ind_j=0;ind_j<signal_length;ind_j++) {
            res_final_real[ind_i+ind_j*window_r2] = K_real_inv_fft_time_1_D[signal_length*ind_i+ind_j];
            res_final_imag[ind_i+ind_j*window_r2] = K_imag_inv_fft_time_1_D[signal_length*ind_i+ind_j];
        }
    }

	// Take fourier transform in Lag domain to form Time-Frequency
	// Distribution using Compact Support Kernel//
    mxSetPr(fft_rhs_final, res_final_real);
    mxSetPi(fft_rhs_final, res_final_imag);
    // Calling Matlab function fft for taking fourier transform in Lag domain
    mexCallMATLAB(1, &fft_lhs_final, 1, &fft_rhs_final, "fft" );
    real_final_time_1_D = mxGetPr(fft_lhs_final);

    for (ind_i=0;ind_i<window_r2*signal_length;ind_i++) {
        result_r [ind_i] = real_final_time_1_D[ind_i];
    }
}


// Function for Calculating Modulus of Indexes
int Modulus_Index(int ind, int N) {
    int out_ind;
    out_ind = ind % N;
    if (out_ind<=0) {
        out_ind=out_ind+N;
    }
    out_ind = out_ind-1;
    return out_ind;
}

// Function for Calculating Analytic Signal
void analytic_signal( double *signal_real,
        double *z_real_out, double *z_imag_out, int signal_length){

    int i;
    double *z_real, *z_imag;
    double *z_real_extend, *z_imag_extend;
    double *z_ifft_real, *z_ifft_imag;
    double *z_t_real, *z_t_imag;

    MATRIX *fft_lhs_data, *rhs_2_data, *fft_rhs_data, *lhs_2_data;

    rhs_2_data = MXCREATEFULL( 1, signal_length*2, COMPLEX );
    fft_rhs_data = MXCREATEFULL( 1, signal_length*2, COMPLEX );

    z_real = (double *) mxMalloc(signal_length*2 * sizeof(double));
    z_imag = (double *) mxMalloc(signal_length*2 * sizeof(double));
    z_real_extend = (double *) mxMalloc(signal_length*2 * sizeof(double));
    z_imag_extend = (double *) mxMalloc(signal_length*2 * sizeof(double));
    z_t_real = (double *) mxMalloc(signal_length*2 * sizeof(double));
    z_t_imag = (double *) mxMalloc(signal_length*2 * sizeof(double));
    z_ifft_real = (double *) mxMalloc(signal_length*2 * sizeof(double));
    z_ifft_imag = (double *) mxMalloc(signal_length*2 * sizeof(double));

    // Extending the length of input signal of length N to 2*N by zero padding
    for (i=0;i<signal_length*2;i++) {
        if (i<signal_length) {
            z_real_extend[i] = signal_real[i];
            z_imag_extend[i] = 0;
        }
        else {
            z_real_extend[i] = 0;
            z_imag_extend[i] = 0;
        }
    }

    mxSetPr(fft_rhs_data, z_real_extend);
    mxSetPi(fft_rhs_data, z_imag_extend);
    // Calling Matlab function fft for computing fourier transform of signal
    mexCallMATLAB(1, &fft_lhs_data, 1, &fft_rhs_data, "fft" );
    z_ifft_real = mxGetPr(fft_lhs_data);
    z_ifft_imag = mxGetPi(fft_lhs_data);

    // Multiplying the first half of the entexded input signal by 2
    for (i=1;i<signal_length;i++) {
        z_ifft_real[i] = 2*z_ifft_real[i];
        z_ifft_imag[i] = 2*z_ifft_imag[i];
    }
    // Replacing the second half of the entexded input signal by 0
    for (i=signal_length+1;i<2*signal_length;i++) {
        z_ifft_real[i] = 0;
        z_ifft_imag[i] = 0;
    }

    for (i=0;i<2*signal_length;i++) {
        z_t_real[i] = z_ifft_real[i];
        z_t_imag[i] = z_ifft_imag[i];
    }
    mxSetPr(rhs_2_data, z_t_real);
    mxSetPi(rhs_2_data, z_t_imag);
    // Calling Matlab function ifft for computing inverse fourier transform
    mexCallMATLAB(1, &lhs_2_data, 1, &rhs_2_data, "ifft" );
    z_real = mxGetPr(lhs_2_data);
    z_imag = mxGetPi(lhs_2_data);

    for (i=0;i<2*signal_length;i++) {
        z_real_out[i] = z_real[i];
        z_imag_out[i] = z_imag[i];
    }
}
