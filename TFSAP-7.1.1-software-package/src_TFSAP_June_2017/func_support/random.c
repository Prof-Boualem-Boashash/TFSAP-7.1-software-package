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
* random.c
*
* random support function
* 
*************************************************/

#include "random.h"

/* define constants as per Table III, L'Ecuyer */

#define	M1 2147483563
#define	M2 2147483399

#define	A1 40014
#define	A2 40692

#define	Q1 53668
#define	Q2 52774

#define	R1 12211
#define	R2 3791

/* state variables */
/* s1, s2 must be in the range [1, M1-1] and [1, M2-1],	respectively */

static long s1 = 123456789;
static long s2 = 398492932;

double rand_gen()
{
  long k, z;

  k = s1 / Q1;
  s1 = A1 * (s1	- k * Q1) - k *	R1;
  if (s1 < 0)
    s1 += M1;

  k = s2 / Q2;
  s2 = A2 * (s2	- k * Q2) - k *	R2;
  if (s2 < 0)
    s2 += M2;

  z = s1 - s2;
  if (z	< 1) z += M1 - 1;

  return (double)( (double)z / (double)M1 );
}


int rand_seed( long init1, long	init2)
{

  if (init1<1)
    init1 = 1;
  if (init2 < 1)
    init2 = 1;
  if (init1>=M1)
    init1 = M1 - 1;
  if (init2>=M2)
    init2 = M2 - 1;

  s1 = init1;
  s2 = init2;

  return 0;
}



