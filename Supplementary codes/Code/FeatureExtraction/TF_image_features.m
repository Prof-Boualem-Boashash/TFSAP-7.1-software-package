function F=TF_image_features(tfd,option)
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
% Time-Frequency Feature Based on Image Processing Techniques
%
% input:
%   tfd --> time-frequency representation
% output:
%   F --> image-related feature vector
%

if option==1
    s  = regionprops(tfd,'Area');
    F(1) = s.Area; % Area or Convex Hull
    
    s  = regionprops(tfd,'Perimeter');
    F(2) = s.Perimeter; % Perimetre
    F(3) = F(2)^2/F(1); % Compactness
    s  = regionprops(tfd,'Centroid');
    F(4) = s.Centroid(2); % Y coordinate of centred region
    F(5) =s.Centroid(1); % X coordinate of centred region
else
    F(1) = image_moments(tfd,0,0); % Area or Convex Hull
    F(2) = sum(sum(bwperim(tfd,8))); % Perimetre
    F(3) = F(2)^2/F(1); % Compactness
    F(4) = image_moments(tfd,1,0)/F(1); % Y coordinate of centred region
    F(5) = image_moments(tfd,0,1)/F(1); % X coordinate of centred region
end

function [mo]=image_moments(tfd,p,q) % Image graysacle moments
[N,M]=size(tfd);
mo=0;

s  = regionprops(tfd,'Centroid');

for i=1:N
    for j=1:M
        mo=mo+((i-s.Centroid(1))^p)*((j-s.Centroid(2))^q)*tfd(i,j);
    end
end