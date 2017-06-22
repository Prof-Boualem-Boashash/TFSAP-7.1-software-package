function Hfeatures = HaralickFeatAllReduced(I, offset, numLevels, grayLimits, compMaxCorrCoef)
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
% Authors: % Boualem Boashash         (boualem.boashash@gmail.com)
%           Samir Ouelha  			(samir_ouelha@hotmail.fr)
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
% HaralickFeatReduced computes the reduced set of Haralick features
% according to Haralick1973 paper
% HARALICKFEATALLREDUCED computes all Haralick features in one step

% Compute GLCMs
offsets = [0 offset;
    -offset offset;
    -offset 0;
    -offset -offset];
glcms = graycomatrix(I, 'Offset', offsets, ...
                        'NumLevels', numLevels, ...
                        'GrayLimits', grayLimits, ...
                        'Symmetric', true);

features = HaralickFeatAll(glcms, compMaxCorrCoef);
Hfeatures = HaralickFeatReduced(features);

end

