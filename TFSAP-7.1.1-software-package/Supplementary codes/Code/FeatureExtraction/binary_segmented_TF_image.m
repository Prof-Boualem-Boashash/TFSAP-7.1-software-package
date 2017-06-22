function [L]=binary_segmented_TF_image(tfd_img, jMax)

%%*************************************************
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

% Segmentation of the Time-Frequency Image
% (Methodoloy based on the watershed algorithm)
%
% input:
%   TFD --> time-frequency representation
% output:
%   bs_img --> binary-segmented T-F image
%


tfd_img=tfd_img/max(max(tfd_img));
T=adaptive_thresh_image(tfd_img,jMax);
tfd_img(tfd_img<T)=0;
tfd_img(tfd_img>=T)=1;
[D,N]=bwlabel(tfd_img);
a=zeros(1,N);

for k=1:1
    for i=1:N     
        a(i)=sum(find(D==i));  
    end
    L=zeros(size(tfd_img));   
    L(D==find(a==max(a)))=1;
    a(a==max(a))=0;
end
end