function ExtractedFeatures=featureComputation(tfd,sig_cur_epoch)
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

% signal related features
cfeatures_vector = signal_featuresExtraction(tfd);

% LBP features
If = myLBP(tfd);
T1=mean(If(:));
T2=std(If(:));
T3=T2/T1;
T4=skewness(If(:));
T5=kurtosis(If(:));
T6= entropy(If);
lbpFeatures = [T1 T2 T3 T4 T5 T6];

%Haralick Features
Hfeatures = HaralickFeatAllReduced(tfd, 1, 16, [], 0);

% Discrete wavelet transform features
dwtFeatures = extractDwtFeatures(sig_cur_epoch,'db4');
waveletFeatures=dwtFeatures ;

ExtractedFeatures = [waveletFeatures lbpFeatures cfeatures_vector Hfeatures];