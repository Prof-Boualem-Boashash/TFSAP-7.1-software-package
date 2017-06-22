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
* gateway routine for cmpt.c
*************************************************/
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mex.h"
#include "matrix.h"
#include "arithm.h"		     /* routines for complex data */
#include "tlocal.h"		     /* local function prototypes */

#include "tfsa_c.h"
#include "cmpt.h"


void mexFunction( int nlhs, MATRIX *plhs[], int nrhs, NEED_CONST MATRIX *prhs[]) {
    int i, x,j;
    double *signal_real, *theta, *E1, *D1;
    int S_M, S_N, theta_M, theta_N, D_M, D_N, E_M,E_N;
    double C,D, E;
    char *p;
    int n;
    int signal_length, theta_length, D_length, E_length, window_r2;
    int window_order;
    int kernel_type;
    double *result_r;
    double **G_real;

    if( nrhs < 4 ) {
        tfsaErr( "Compacttfd", "Not enough input arguments" );
        winNTcheck( nlhs, plhs );
        return;
    }

    if( nlhs > 1 ) {
        tfsaErr( "Compacttfd", "Too many output arguments" );
        winNTcheck( nlhs, plhs );
        return;
    }



    S_M = (int)mxGetM( prhs[0] );
    S_N = (int)mxGetN( prhs[0] );
    if( (S_M==1 && S_N==1) || (S_M!=1 && S_N!=1) ) {
        mexPrintf("\n Invalid Dimensions ");
    }
    /* Choose the larger row or column dimension for signal length */
    signal_length = (unsigned)( S_M>S_N ? S_M:S_N );

    // Real Input Signal
    signal_real = mxGetPr(prhs[0]);

    window_order = 0;
    window_r2 = 1;
    while( window_r2 < signal_length) {
        window_order++;
        window_r2 <<= 1;
    }

    result_r = (double *) mxMalloc( window_r2 * signal_length * sizeof(double));


    /* Fourth input argument */

    if( MXSTRING( prhs[1] ) ) {

        /* Get string from Matlab following p2-51 in Mex manual */

        n = (int)mxGetN( prhs[1] ) + 1;  /* string length + 1 for NULL */
        p = mxCalloc( n, sizeof( char ) );
        if( p == NULL ) {
            tfsaErr( "Compacttfd", "Internal memory allocation failure" );
            winNTcheck( nlhs, plhs );
            return;
        }

        if( mxGetString( prhs[1], p, n ) ) {
            tfsaErr( "Compacttfd", "Could not get kernel type" );
            winNTcheck( nlhs, plhs );
            return;
        }

        if( !strcmp( p, "csk" ) || !strcmp( p, "CSK" ) )
            kernel_type = CSK;
        else if( !strcmp( p, "ckd" ) || !strcmp( p, "CKD" ) )
            kernel_type = CKD;
        else {
            tfsaErr( "Compacttfd", "Unknown kernel type" );
            winNTcheck( nlhs, plhs );
            return;
        }
    }
    else {
        tfsaErr( "Compacttfd", "Smoothing window type must be a string" );
        winNTcheck( nlhs, plhs );
        return;
    }


    G_real = (double **)mxCalloc( window_r2, sizeof(double *) );
    if( G_real == NULL ) {
        mexPrintf( "\nMemory allocation	failed\n"	);
    }
    G_real[0] = (double *) mxCalloc( (int)(window_r2*signal_length), sizeof(double));
    if( G_real[0] == NULL) {
        mexPrintf( "\nMemory allocation	failed\n"	);
    }
    /* Set up pointers for pointer array */
    for(i = 1; i < window_r2; i++ ) {
        G_real[i] = G_real[0] + (unsigned long)signal_length * i;
    }
    for (i=0;i<window_r2;i++) {
        for (j=0;j<signal_length;j++) {
            G_real[i][j] = 0;
        }
    }

    
    switch( kernel_type ) {

        case( CSK ):

            //// ------ % Filter Design Parameters for Compact Support Kernel ----//////
            // Filter Parameters C, D
            
			for (x=2;x<=3;x++) {
               
				if (mxIsChar(prhs[x]) || !mxIsDouble(prhs[x]) || mxIsComplex(prhs[x]) ||
               mxGetN(prhs[x])*mxGetM(prhs[x]) != 1) 
				{
                  mexPrintf("Inputs must be a scalar.");
			      return;
                }
	}
			C = *mxGetPr(prhs[2]);
            D = *mxGetPr(prhs[3]);


            csk_kernel(G_real, window_r2, signal_length, C, D);

            break;

        case( CKD ):

            if( nrhs < 5 ) {
                tfsaErr( "Compacttfd", "Not enough input arguments" );
                winNTcheck( nlhs, plhs );
                return;
            }

            //// ------ % Filter Design Parameters for Extended Compact Support Kernel ----//////
            // Filter Parameters C, D, E
            
			for (x=2;x<=4;x++) {
               
				if (mxIsChar(prhs[x]) || !mxIsDouble(prhs[x]) || mxIsComplex(prhs[x]) ||
               mxGetN(prhs[x])*mxGetM(prhs[x]) != 1) 
				{
                  mexErrMsgTxt("Inputs must be a scalar.");
			      return;
                }
	        }
			
			C = *mxGetPr(prhs[2]);
            D = *mxGetPr(prhs[3]);
            E = *mxGetPr(prhs[4]);
			


			  
            ecsk_kernel(G_real, window_r2, signal_length, C, D, E);//

            break;
    }

    // Calling cmpt_filter function for performing filtering
    cmpt_filtering( signal_real, signal_length,
            window_r2, G_real, result_r);

    plhs[0] = MXCREATEFULL(  window_r2, signal_length, REAL );

    if (plhs[0] == NULL) {
        mexPrintf("Memory allocation failure: plhs");
    }

    mxSetPr(plhs[0], result_r);

}


void csk_kernel(double **G_real, int Rows, int Columns, double C, double D) {
    int ii, jj,offset;
    double val, val_i, val_j;

    for (ii=0;ii<Rows/2;ii++) {
        for (jj=0;jj<floor(Columns/2);jj++) {
            val_i = ((double)(ii))*1/Rows;
            val_j = ((double)(jj))*1/Columns;

            if ((pow(val_i, 2)<pow(D, 2)) && (pow(val_j, 2)<pow(D, 2)) ) {
            if (abs(val_i-val_j)<D)
                val = exp(2*C)*exp(C*pow(D, 2)*(1/(pow(val_i, 2)-pow(D, 2))+1/(pow(val_j, 2)-pow(D, 2))));

                G_real[ii][jj]= val;
                G_real[Rows - 1- ii][jj] = val;
                if (jj>0) {
                    G_real[ii] [Columns  - jj] = val;
                    G_real[Rows - 1 - ii] [Columns  - jj] = val;
                }
            }
        }
    }



}

void ecsk_kernel(double **G_real, int Rows, int Columns, double C, double D, double E) {

    int ii, jj, offset; 
    double val, val_i, val_j, chk_val;
	//FILE *fichier = fopen("check.txt", "w" );

    for (ii=0;ii<Rows/2;ii++) {
        for (jj=0;jj<floor(Columns/2);jj++) {
            val_i = ((double)(ii))*1/Rows;
            val_j = ((double)(jj))*1/Columns;

            chk_val = pow(val_i, 2)/pow(D, 2)+ pow(val_j, 2)/pow(E, 2);

            if (chk_val<1) {
                val = exp(2*C)*exp((C*pow(D, 2)/(pow(val_i, 2)-pow(D, 2)))+C*pow(E, 2)/(pow(val_j, 2)-pow(E, 2)));
                G_real[ii][jj]= val;
                G_real[Rows - 1- ii][jj] = val;
                if (jj>0) {
                    G_real[ii] [Columns  - jj] = val;
                    G_real[Rows - 1 - ii] [Columns  - jj] = val;
                }
            }
        }
    }

}



 

                


