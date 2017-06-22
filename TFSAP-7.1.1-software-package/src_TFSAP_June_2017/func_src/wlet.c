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
* Wavelet Transform
* Performs the	forward	or inverse wavelet transform using a pyramid algorithm.
* The transform is performed in-place.
*************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "wlet.h"

/* Daubechies wavelet coefficients */

double c4[5] = {
		   0.4829629131445341, 0.8365163037378079,
		   0.2241438680420134, -0.1294095225512604};

double c4r[5] =	{  -0.1294095225512604,	-0.2241438680420134,
		   0.8365163037378079, -0.4829629131445341};

double c12[13] = {
		   0.111540743350, 0.494623890398, 0.751133908021,
		   0.315250351709,-0.226264693965,-0.129766867567,
		   0.097501605587, 0.027522865530,-0.031582039318,
		   0.000553842201, 0.004777257511,-0.001077301085};

double c12r[13]	= {-0.001077301085,-0.004777257511, 0.000553842201,
		   0.031582039318, 0.027522865530,-0.097501605587,
		  -0.129766867567, 0.226264693965, 0.315250351709,
		  -0.751133908021, 0.494623890398,-0.111540743350 };

double c20[21] = {
		   0.026670057901, 0.188176800078, 0.527201188932,
		   0.688459039454, 0.281172343661,-0.249846424327,
		  -0.195946274377, 0.127369340336, 0.093057364604,
		  -0.071394147166,-0.029457536822, 0.033212674059,
		   0.003606553567,-0.010733175483, 0.001395351747,
		   0.001992405295,-0.000685856695,-0.000116466855,
		   0.000093588670,-0.000013264203};

double c20r[21]	= { -0.000013264203,-0.000093588670,-0.000116466855,
		    0.000685856695, 0.001992405295,-0.001395351747,
		   -0.010733175483,-0.003606553567, 0.033212674059,
		    0.029457536822,-0.071394147166,-0.093057364604,
		    0.127369340336, 0.195946274377,-0.249846424327,
		   -0.281172343661, 0.688459039454,-0.527201188932,
		    0.188176800078,-0.026670057901};

unsigned  n_coeff;
double *coeff_f, *coeff_r;
int off, off;
int t_dir;			/* transform direction */
int r_start;			/* number of points in the stage at
				 which to start	the reverse
				 decomposition,	equal to the smallest
				 power of two greater than the number
				 of coefficients */

/* returns -1 if error,	otherwise the number of	stages.
 */

int wave( double *data,	double *result,	unsigned long n, int direction,	unsigned num_coeff)
{
  unsigned long	step_n;
  unsigned i, stages;
  double *temp;

  if (n	< num_coeff) return -1;

  init_wave( num_coeff,	direction);

  for(i	= 0; i < n; i++)
    result[i] =	data[i];

  temp = (double *)malloc(sizeof(double) * (int)n);
  if (temp == NULL)
    return -1;

  stages = 0;

  if (direction	== 1)
    /* forward transform */
    /* work from largest to smallest */
    for	( step_n = n; step_n >=	n_coeff; step_n	>>= 1) {
      wave_step( result, temp, step_n);
      stages++;
    }
  else
    /* reverse transform */
    /* work from smallest to largest */
    for	( step_n = r_start; step_n <= n; step_n	<<= 1){
      wave_step( result, temp, step_n);
      stages++;
    }

  free(temp);

  return stages;

}

int init_wave( int num_coeff, int direction)
{
  n_coeff = num_coeff;
  t_dir	= direction;

  if (num_coeff	== 20) {
    coeff_f = c20;
    coeff_r = c20r;
    r_start = 32;
  }
  else if (num_coeff ==	12) {
    coeff_f = c12;
    coeff_r = c12r;
    r_start = 16;
  }
  else {
    coeff_f = c4;
    coeff_r = c4r;
    n_coeff = 4;
    r_start = 4;
  }

  /* offsets for the convolution, to centre the	output */
  off =	-((int)n_coeff >> 1)+1;

  return 0;
}


int wave_step( double *data, double *temp, unsigned long n)
{
  unsigned long	i, ii, j, index, k, nD2;
  double d1, d2;

  if (n	< n_coeff)
    return -1;

  nD2 =	n >> 1;

  for (j = 0; j	< n; j++)
    temp[j] = 0.0;

  if (t_dir == 1) {
    for	(ii = 0, i = 0;	i < n; i+= 2, ii++) {
      /* i is the index	into data at which the convolution will	start */
      for (k = 0; k < n_coeff; k++) {
	index =	(i + k + off) %	n;
	temp[ii] += coeff_f[k] * data[index]; /* low pass */
	temp[ii+nD2] +=	coeff_r[k] * data[index]; /* high pass */
      }
    }
  }
  else {
    for	(ii = 0, i = 0;	i < n; i+= 2, ii++) {

      d1 = data[ii];
      d2 = data[ii+nD2];

      for (k = 0; k < n_coeff; k++) {
	index =	(i + k + off) %	n;
	temp[index] += coeff_f[k] * d1;	/* low pass */
	temp[index] += coeff_r[k] * d2;	/* high	pass */
      }
    }
  }
  for (i=0; i<n; i++)
    data[i] = temp[i];

  return 0;
}


/*
 * This	routine	converts 1D array of wavelet coefficients
 * into	the 2D time-frequency array, according to the
 * famous wavelet tiling. Note that the	output is an
 * energy distribution (scalogram).
 */


int form_ts( double *input, double *output, int	n )
{
 int p,	r, N;
 int i,	 k;
 int xind, step, substep, ind1,	ind2, bl, bs, be;
 int two2step, two2step1;

   N = n;
  r = 0;
  p = 1;
  while	(p < N)	{
	r++;
	p = p*2;
	}

/* smoothing coefficients  */
  for (i = 0; i	< N/2; i++) {output[i] = input[0]*input[0];}
  for (i = N/2;	i < N; i++) {output[i] = input[1]*input[1];}

/* wavelet coefficients	*/
  xind = 1;
  for (step=1; step <= r-1; step++) {
	two2step = 1;
	for (i=1; i<=step; i++)	{two2step = two2step * 2;}
	two2step1 = 1;
	for (i=1; i<=step+1; i++) {two2step1 = two2step	* 2;}
	for (substep = 0; substep < two2step; substep ++) {
		ind1 = two2step	+ substep;
		ind2 = N / two2step1;
		for (bl	= 1; bl	<= two2step; bl++) {
			bs = (bl-1)*ind2;
			be = (bl * ind2)-1;
			for (k=bs; k <=	be; k++)
			    {output[N/2*ind1+k]	= input[xind+bl]*input[xind+bl];}

			}
		}
	xind = xind+two2step;
	}
return 0;
}

