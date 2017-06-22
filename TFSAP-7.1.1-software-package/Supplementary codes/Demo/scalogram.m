function [TFD,CWT,f] = scalogram(x,wave, M)
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
%
% Last Modification: 25-12-2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
%
% Computes the continuous wavelet transform (CWT) and the scalogram using
% Morlet wavelt.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialisation and parameters
thr  = 1e-2;
fmin = thr;
fmax = 0.5;
N = length(x);
f = logspace(log10(fmin),log10(fmax),M);
a = logspace(log10(fmax/fmin),log10(1),M);

%% Wavelets generation
s = (real(x) - mean(real(x)))';
z = hilbert(s);
CWT = zeros(M, N);
TFD = zeros(M, N);
if(wave >= 1)
    for k = 1:M
        nha = wave*a(k);
        th  = -round(nha) : round(nha);
        h   = exp(-(2*log(10)/nha^2)*th.^2).*exp(1i*2*pi*f(k)*th); % Morlet Wavelet
        CWT(k,:) = conv(z, h,'same')./sqrt(a(k));
        TFD(k,:) = CWT(k,:).*conj(CWT(k,:)) ;
    end

else
    fprintf(2,'Smoothing parameters should be more than 1 ...\n'); return;
end
