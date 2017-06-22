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
*						  RAs: B.R.  and B.H.
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
*
* Description:  
 *  COMPLEX ARITHMETICS
 *
 *   complex cmplx (double a, double b)   - create a complex number
 *   complex add (complex x, complex y)   - add two complex numbers
 *   complex sub (complex x, complex y)   - substract two complex numbers
 *   complex multpl(complex x, complex y)  - multiply two complex numbers
 *   complex conj (complex x)       - complex-conjugate of complex number
 *   double  modul (complex x)      - absolute value of complex number
 *   double  arg (complex x)        - argument of complex number
 *   complex cmult (complex x, double c)  - muliply by real constant
 *
 *
 *****************************************************************************/

#ifndef ARITHM_H
#define ARITHM_H

typedef struct
       {
       double re;
       double im;
       }           complex;

complex cmplx (double a, double b);
complex add (complex x, complex y);
complex sub (complex x, complex y);
complex multpl(complex x, complex y);
complex divide(complex x, complex y);
complex conj (complex x);
double  module (complex x);
double  arg (complex x);
complex cmult (complex x, double c);

#endif

