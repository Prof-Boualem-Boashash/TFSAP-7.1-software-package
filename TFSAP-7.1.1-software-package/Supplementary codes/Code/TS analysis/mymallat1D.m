function[coeffs,dimplan]=mymallat1D(h,g,signal,dimplan)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 Boualem Boashash
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%
% Authors: Prof. Boualem Boashash   (boualem.boashash@gmail.com)
%          Samir Ouelha             (samir_ouelha@hotmail.fr)
%          
% The following references should be cited whenever this script is used:
% [1] B. Boashash, Samir Ouelha, Designing high-resolution time-frequency
% and time-scale distributions for the analysis and  classification of non 
% stationary signals: a tutorial review with features performance comparison
%, Digital Signal Processing, In Press.
% [2] B. Boashash, Samir Ouelha, Efficient Software Matlab Package to compute
% Time-Frequency Distributions and related Time-Scale methods with extraction
% of signal characteristics, SoftwareX, In Press.
%
% This study was funded by grants from the ARC and QNRF NPRP 6-885-2-364
% and NPRP 4-1303-2-517 
%% Last Modification: 25-12-2016
%
% Description:
%
%This function uses the pyramidal multiresolution algorithm to compute DWT
% coefficients.
%
% Inputs: h = Impulse response of the low pass filter.
% g = impulse response of the high pass filter
% signal = input signal if nargin = 3 or coefficients if nargin 4 
% dimplan: dimension of each level of decomposition for the synthesis of
% the signal.
%
%Outputs: 1) if nargin=3, coeff represents the DWT coefficients and if
%nargin=4 coeff represents the synthesis of the input coefficients.
% dimplan: represents the dimension of each level of decomposition.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

signal = signal(:)';

if nargin==3
    N = length(signal);
    Nh = length(h);
    res = ceil(log2(N/Nh));
    
    htild = h(end:-1:1);
    gtild = g(end:-1:1);
    
    coeffs = [ ];
    dimplan = zeros(1,res+2);
    
    approx = signal;
    
    for n = 1 : res
        a = myconv(approx,htild);
        d = myconv(approx,gtild);
        a = a(2:2:end);
        d = d(2:2:end);
        coeffs = [coeffs,d]; %#ok<AGROW>
        approx = a;
        dimplan(n) = length(d);
    end
    
    coeffs = [coeffs,a];
    dimplan(res+1) = length(a);
    dimplan(res+2) = N;
    
else 
    Nh = length(h);
    res = length(dimplan)-2;
    N = dimplan(end);
    dimplan = dimplan(1:end-1);
    a = signal(end-dimplan(end)+1:end);
    fin = dimplan(end);
    
    for n = 1:res
        d = signal(end-fin-dimplan(end-n)+1:end-fin);
        fin = fin + dimplan(end-n);
        
        areech = zeros(1,2*length(a));
        dreech = zeros(1,2*length(d));
        
        areech(1:2:end) = a;
        dreech(1:2:end) = d;
        
        areech = myconv(areech(1:end-1),h); % one zero is removed
        dreech = myconv(dreech(1:end-1),g); % one zero is removed
        
        a = areech + dreech;
        
        if n < res
            a = a(Nh-1:Nh+dimplan(end-n-1)-2);
        else 
            a = a(Nh-1:Nh+N-2);
        end
    end
    
    coeffs = a;
    
end