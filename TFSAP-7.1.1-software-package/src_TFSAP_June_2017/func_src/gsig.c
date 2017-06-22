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
* Generation of various test signals
*************************************************/

#include <stdlib.h>
#include <math.h>
#include "arithm.h"
#include "fft.h"
#include "gsig.h"


#define	TFSA_PI	3.14159265358979323846




#define	EPS 1e-5			/* arbitrary, but small	offset */

#ifdef DLLMEX							/* Add offset for PC version only... */
#define	COS(x)	cos( (double)(x	+ EPS) )
#define	SIN(x)	sin( (double)(x	+ EPS) )
#else
#define	COS(x) cos(x)
#define	SIN(x) sin(x)
#endif




void linfm (double *result_real, double	*result_imag, double f1, double	f2, int	n)
{
  int i;
  double phi;

  if( result_imag == NULL )
    for	(i = 0;	i < n; i++) {
      phi = (f1	* (double)i + 0.5 * (f2-f1) * (double)i	* (double)i / (double)n	);

      result_real[i] =	  COS( 2.0*TFSA_PI*phi );
    }
  else
    for	(i=0; i<n; i++)	{
      phi = ( f1*(double)i + (f2-f1)*0.5*(double)i*(double)i/(double)n );
      result_real[i] = COS( 2.0*TFSA_PI*phi );
      result_imag[i] = SIN( 2.0*TFSA_PI*phi );
    }
}


void sinfm( double *result_real, double	*result_imag, double fc,
	   double fm, int n, double fdev )
{
  int i;
  double phi;	     /*	phase */

  if (result_imag == NULL)
    for	(i=0; i<n; i++)	{
      phi = (fc*(double)i) + fdev/2.0*SIN( 2.0*TFSA_PI*fm*(double)i );
      result_real[i] =	COS( 2.0*TFSA_PI*phi );
    }
  else
    for	(i=0; i<n; i++)	{
      phi = (fc*(double)i)+ fdev/2.0*SIN( 2.0*TFSA_PI*fm*(double)i );
      result_real[i] = COS( 2.0*TFSA_PI*phi );
      result_imag[i] = SIN( 2.0*TFSA_PI*phi );
    }
}


void quadfm (double *result_real, double *result_imag, double f1, double f2, int n)
{
  int i;
  double phi;

  if (result_imag == NULL)
    for	(i=0; i<n; i++)	{
      phi = ( f1*(double)i + 2.*(f2-f1)*(double)i*(double)i/n -	4.*(f2-f1)*(double)i*(double)i*(double)i/(3.*n*n) );
      result_real[i] = COS( 2*TFSA_PI*phi );
       }
  else
    for	(i=0; i<n; i++)	{
      phi = ( f1*(double)i + 2.*(f2-f1)*(double)i*(double)i/n -	4.*(f2-f1)*(double)i*(double)i*(double)i/(3.*n*n) );
      result_real[i] = COS( 2.0*TFSA_PI*phi );
      result_imag[i] = SIN( 2.0*TFSA_PI*phi );
    }
}



void hypfm( double *result_real, double	*result_imag, double f1,
	   double f2, int n)
{
  int i;
  double phi;
  int type = 1;

  double beta0,	f0;  /*	parameters that	correspond to central and
		       "scale" parameters */

  /*
   * Knowing F1	(starting frequency) and F2 (end frequency)
   * you can calculate f0 and beta0 parameters imposing
   * equations:	F1 = f0/(1-beta0*n/2); F2 = f0/(1+beta0*n/2)
   * The inst. frequency is: fi	= f0/(1.+beta0*(double)(i-n/2))
   *
   */

  beta0	= (2.*(f1-f2)) / ((f1+f2)*(double)n);
  f0 = f1 * (1.	- beta0	* (double)n / 2.);

  if( result_imag == NULL)
    for	(i= 0; i<n; i++)
      {
	if (type == 1)
	  {
	    phi	= (f0/beta0)*log(fabs(1.+ beta0	* ((double)i-(double)n/2.0)));
	    result_real[i] = COS( 2.0*TFSA_PI*phi );
	  }
	else
	  {
	    phi	= f1 * log ( (double)i+1.);
	    result_real[i] = COS( 2.0*TFSA_PI*phi );
	  }
      }
  else
    for	(i= 0; i<n; i++)
      {
	if (type == 1)
	  {
	    phi	= (f0/beta0)*log(fabs(1.+ beta0	* ((double)i-(double)n/2.0)));
	    result_real[i] = COS( 2.0*TFSA_PI*phi );
	    result_imag[i] = SIN( 2.0*TFSA_PI*phi );
	  }
	else
	  {
	    phi	= f1 * log ( (double)i+1.0);
	    result_real[i] = COS( 2.0*TFSA_PI*phi );
	    result_imag[i] = SIN( 2.0*TFSA_PI*phi );
	  }
      }
}



void cubicfm( double *result_real, double *result_imag,	double f0, double f1, int n)
{
  int i;
  double phi;
  double cf = 0.33333333333;	/* centre frequency (where max occurs) */
  double a, b, c, d;

  /* the equation for phi was determined by solving the	following constraints:
     >>	[x,y,z]	= solve('a*cf^3*n^3 + b*cf^2*n^2 + c*cf*n + f0 = (f1-f0)/2 + f0', ...
     'a	* n^3 +	b * n^2	+ c * n	+ f0 = f1', '3*a*cf^2*n^2+2*b*cf*n+c=0', 'a,b,c');

     x =

     -1/2/cf^2*(2*f1*cf+2*cf^2*f0-2*cf^2*f1-2*f0*cf+f0-f1)/n^3/(-2*cf+cf^2+1)

     >>	y

     y =

     1/2*(3*cf^2*f1+f0+4*cf^3*f0-4*cf^3*f1-f1-3*cf^2*f0)/cf^2/n^2/(-2*cf+cf^2+1)

     >>	z

     z =

     -1/2*(2*cf^3*f0-2*cf^3*f1-3*f0*cf+3*f1*cf+2*f0-2*f1)/cf/n/(-2*cf+cf^2+1)

     */

  a = -1.0/2.0/cf/cf*(2.0*f1*cf+2.0*cf*cf*f0-2.0*cf*cf*f1-2.0*f0*cf+f0-f1)/(double)n/(double)n/(double)n/(-2.0*cf+cf*cf+1.0);
  b = 1.0/2.0*(3.0*cf*cf*f1+f0+4.0*cf*cf*cf*f0-4.0*cf*cf*cf*f1-f1-3.0*cf*cf*f0)/cf/cf/(double)n/(double)n/(-2.0*cf+cf*cf+1.0);
  c = -1.0/2.0*(2.0*cf*cf*cf*f0-2.0*cf*cf*cf*f1-3.0*f0*cf+3.0*f1*cf+2.0*f0-2.0*f1)/cf/(double)n/(-2.0*cf+cf*cf+1.0);
  d = f0;

  if( result_imag == NULL)
    for	(i=0; i<n; i++)	{
      phi = a *	(double)i * (double)i *	(double)i * (double)i /	4.0;
      phi += b * (double)i * (double)i * (double)i / 3.0;
      phi += c * (double)i * (double)i / 2.0;
      phi += d * (double)i;
      result_real[i] = COS( 2.0*TFSA_PI*phi );
    }
  else
    for	(i=0; i<n; i++)	{
      phi = a *	(double)i * (double)i *	(double)i * (double)i /	4.0;
      phi += b * (double)i * (double)i * (double)i / 3.0;
      phi += c * (double)i * (double)i / 2.0;
      phi += d * (double)i;
      result_real[i] = COS( 2.0*TFSA_PI*phi );
      result_imag[i] = SIN( 2.0*TFSA_PI*phi );
    }
}


void stepfm( double *result_real, double *result_imag, double f1, double f2, int n,
	    int	num_bins)
{
  int i;
  double f;
  double increment = (double)(f2-f1) / (num_bins -1);
  double bin_size = (double)n /	(double)num_bins;


  double sum = 0.0;
  double wind[WINLEN] =	{ 0, 0,	0};
  double smoothed_freq;


  if (result_imag == NULL) {
    for	(i = 0;	i < num_bins * bin_size; i++) {
      f	= f1 + (int)(i / bin_size) * increment;

      if (i < WINLEN) {
	smoothed_freq =	f;
	sum += f;
	wind[i]	= f;
      }
      else {
	sum += f - wind[i%WINLEN];
	smoothed_freq =	sum/WINLEN;
	wind[i%WINLEN] = f;
      }
      result_real[i] = COS( 2.0*TFSA_PI*smoothed_freq*(double)i	);
    }

    for( ; i<n;	i++) {
      result_real[i] = COS( 2.0*TFSA_PI*smoothed_freq*(double)i	);
    }
  }
  else {
    for	(i = 0;	i < num_bins * bin_size; i++) {
      f	= f1 + (int)(i / bin_size) * increment;

      if (i < WINLEN) {
	smoothed_freq =	f;
	sum += f;
	wind[i]	= f;
      }
      else {
	sum += f - wind[i%WINLEN];
	smoothed_freq =	sum/WINLEN;
	wind[i%WINLEN] = f;
      }
      result_real[i] = COS( 2.0*TFSA_PI*smoothed_freq*(double)i	);
      result_imag[i] = SIN( 2.0*TFSA_PI*smoothed_freq*(double)i	);
    }

    for( ; i<n;	i++) {
      result_real[i] = COS( 2.0*TFSA_PI*smoothed_freq*(double)i	);
      result_imag[i] = SIN( 2.0*TFSA_PI*smoothed_freq*(double)i	);
    }
  }


}

