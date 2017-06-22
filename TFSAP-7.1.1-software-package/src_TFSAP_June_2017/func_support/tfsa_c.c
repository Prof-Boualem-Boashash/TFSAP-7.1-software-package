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
* tfsa_c.c: support functions for TFSAP
* 
*************************************************/
#include <stdio.h>
#include "mex.h"
#include "tfsa_c.h"


/* Global variables for	waitbar*/
static MATRIX *prhs[1];
static MATRIX *p1;
static MATRIX *plhs[1];



void tfsa_cerr(	char *error_txt)
{

  p1 = mxCreateString(error_txt);
  prhs[0] = p1;
  mexCallMATLAB( 0, plhs, 1, prhs, "tfsa_err");

  return;
}

void tfsa_cwarn( char *warn_txt)
{

  p1 = mxCreateString(warn_txt);
  prhs[0] = p1;
  mexCallMATLAB( 0, plhs, 1, prhs, "tfsa_wrn");

  return;
}


double * Win_NT_Err (int nlhs, MATRIX *plhs[])

/*
 * This	is a hack to overcome an error in Matlab mex calls under Windows NT.
 * In windows NT if a functions	is called with one or more lhs arguements
 * and a TFSA error occurs then	the function will cause	a general protection fault.
 * For example	"[x,y]=wvd" causes a TFSA error	because	there is not enough inputs.
 * When	the function is	returned a GPF occurrs.	 To stop this I	have allocated
 * a matrix of size one	for each lhs arguement and intialised it to zero.
 *
 * This	function should	be called at the beginning of each interface mex function.
 *
 */
{
   double * result;
   int i;

   if (nlhs > 0) {
      for (i=0;	i < nlhs; i++) {
	 plhs[i] = MXCREATEFULL( 1, 1, REAL  );
	 result	= mxGetPr (plhs[i]);
	 result[i] = 0;
      }
   }
   return result;
}
