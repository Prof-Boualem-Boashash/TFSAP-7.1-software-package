function dwtFeatures = extractDwtFeatures(signal,wavname)
%*************************************************
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
% Authors: Boualem Boashash         (boualem.boashash@gmail.com)
%          Samir Ouelha  			(samir_ouelha@hotmail.fr)
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
%
% Description:
%
%This function computes statistical features from the DWT decomposition.
%Inputs: 1) signal: this is the signal from where we extract features
%        2)  wavname: represents the wavelet used for the analysis.
%Output: 1) dwtFeatures: represents the statistical features from the DWT.


[LO_R,HI_R] = mywfilters(wavname);

[signal_dec, dimplan] = mymallat1D(LO_R,HI_R,signal);
dwtFeatures=[];
deb = 1;
for n=1:length(dimplan)-2,
    signal_dim=signal_dec(deb:deb+dimplan(n)-1);
    deb = deb+dimplan(n);
    T1=mean(signal_dim);
    T2=std(signal_dim);
    T3=T2/T1;
    T4=skewness(signal_dim);
    T5=kurtosis(signal_dim);
    T6=iqr(signal_dim);
    T = [T1 T2 T3 T4 T5 T6];
    dwtFeatures = [ T dwtFeatures];
end