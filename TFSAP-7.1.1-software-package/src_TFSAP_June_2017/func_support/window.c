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
* window
*
* Windowing functions for use with TFSA
* 
*************************************************/

#include <math.h>
#define	TFSA_PI	3.14159265358979323846



#define	EPS 1e-5			/* arbitrary, but small	offset */

#ifdef DLLMEX							/* Add offset for PC version only... */
#define	COS(x)	cos( (double)(x	+ EPS) )
#define	SIN(x)	sin( (double)(x	+ EPS) )
#else
#define	COS(x) cos(x)
#define	SIN(x) sin(x)
#endif




void bartlett( double *array, int length)
{
  int i;
  float	centre;

  if (length ==	1)
    *array = 1.0;
  else {
    centre = (float)(length-1)/2;
    for	(i = 0;	i< length; i++)
      array[i] = 1 - fabs((i - centre)/centre);
  }
}

void hann( double *array, int length)
{
  int i;

  if (length ==	1)
    *array = 1.0;
  else
    for	(i = 0;	i< length; i++)
      array[i] = 0.5 * (1 - COS( 2.0 * TFSA_PI * (double)i / (double)(length-1)));

}

void hamming( double *array, int length)
{
  int i;

  if (length ==1)
    *array = 1.0;
  else
    for	(i = 0;	i< length; i++)
      array[i] =  (0.54	- 0.46 * COS( 2.0 * TFSA_PI * (double)i	/ (double)(length-1)));

}
