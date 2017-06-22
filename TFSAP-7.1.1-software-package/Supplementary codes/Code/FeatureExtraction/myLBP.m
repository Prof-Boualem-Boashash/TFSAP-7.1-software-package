function I_fin = myLBP(I)
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
% This function transforms an image applying the operator LBP. For more
% information about this paper, please see the paper associated to this
% code.
%
%Input: The image to transform
%Output: The resulting image

I = [zeros(1,size(I,2));I;zeros(1,size(I,2))];
I = [zeros(size(I,1),1),I,zeros(size(I,1),1)];
for nb_i= 2:size(I,1)-1
    for nb_j = 2:size(I,2)-1
        threshold = I(nb_i,nb_j);
        tab=zeros(3,3);
        if I(nb_i-1,nb_j-1)>=threshold
            tab(1,1) = 1;
        end
        if I(nb_i-1,nb_j)>=threshold
            tab(1,2) = 1;
        end
        if I(nb_i-1,nb_j+1)>=threshold
            tab(1,3) = 1;
        end
        if I(nb_i,nb_j+1)>=threshold
            tab(2,3) = 1;
        end
        if I(nb_i+1,nb_j+1)>=threshold
            tab(3,3) = 1;
        end
        if I(nb_i+1,nb_j)>=threshold
            tab(3,2) = 1;
        end
        if I(nb_i+1,nb_j-1)>=threshold
            tab(3,1) = 1;
        end
        if I(nb_i,nb_j-1)>=threshold
            tab(2,1) = 1;
        end
        I(nb_i,nb_j) = tab(1,3)*2^7+tab(1,2)*2^6+tab(1,3)*2^5+tab(2,3)*2^4+tab(3,3)*2^3+tab(3,2)*4+tab(3,1)*2+tab(2,1);
    end
end
I_fin = I(2:end-1,2:end-1);