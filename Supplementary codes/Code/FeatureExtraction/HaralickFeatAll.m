function features = HaralickFeatAll(glcms_in, compMaxCorrCoef)
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
% Authors:  Boualem Boashash         (boualem.boashash@gmail.com)
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
% HARALICKFEATALL computes a set of features proposed by Haralick and others
% Input: glcms : a 3D matrix of 4 glcms
% Output
% Features: The set of Haralick features

if (~isfloat(glcms_in))
    error ('ERROR. glcms must be a floating point matrix');
end

if (nargin == 1)
    compMaxCorrCoef = 0; % To avoid numerical instability
    % of maximal correlation coefficient
elseif (nargin == 2)
    if (~isscalar(compMaxCorrCoef))
        error('ERROR. compMaxCorrCoef should be a scalar 1 or 0');
    end
end
    


glcms = glcms_in;

L = size(glcms(:,:,1), 1);
% glcm
Nglcms = size(glcms, 3);

% Auxiliary variables
glcmSum  = zeros(Nglcms,1);
glcmMean = zeros(Nglcms,1);
glcmVar  = zeros(Nglcms,1);

% Marginal-probabilities
p_x = zeros(L,Nglcms);  
p_y = zeros(L,Nglcms);
p_xplusy = zeros((2*L - 1),Nglcms);
p_xminusy = zeros((L),Nglcms);

mu_x = zeros(Nglcms,1);
mu_y = zeros(Nglcms,1);
sigma_x = zeros(Nglcms,1);
sigma_y = zeros(Nglcms,1);
% Just to double check (another way to do the computations)
mu_x2 = zeros(Nglcms, 1);
mu_y2 = zeros(Nglcms, 1);
sigma_x2 = zeros(Nglcms, 1);
sigma_y2 = zeros(Nglcms, 1);

% Information Measures of Correlation auxiliary variables
hx   = zeros(Nglcms,1); % Entropy of p_x
hy   = zeros(Nglcms,1); % Entropy of p_y
hxy  = zeros(Nglcms,1); % Same as the entropy
hxy_1 = zeros(Nglcms,1);
hxy_2 = zeros(Nglcms,1);

% Maximal Correlation Coefficient auxiliary variables
if (compMaxCorrCoef)
    Q = zeros(size(glcms));
end

% normalization/mean/variance for glcms
for i=1:Nglcms
    glcmSum(i) = sum(sum(glcms(:,:,i)));
    glcms(:,:,i) = glcms(:,:,i) / glcmSum(i);
    glcmMean(i) = mean2(glcms(:,:,i));
    glcmVar(i) = (std2(glcms(:,:,i)))^2;
end

% compute marginal probabilities and their mean/standard deviation
p_xplusy_range = 2:2*L;
p_xminusy_range = 0:(L - 1);
for k = 1:Nglcms
    for i = 1:L
        for j = 1:L
            p_x(i,k) = p_x(i,k) + glcms(i,j,k); 
            p_y(i,k) = p_y(i,k) + glcms(j,i,k);
            if (ismember((i + j),p_xplusy_range))
                p_xplusy((i+j)-1,k) = p_xplusy((i+j)-1,k) + glcms(i,j,k);
            end
            if (ismember(abs(i-j),p_xminusy_range))
                
                p_xminusy((abs(i-j))+1,k) = p_xminusy((abs(i-j))+1,k) ...
                + glcms(i,j,k);
            end
            mu_x(k)          = mu_x(k) + (i)*glcms(i,j,k);
            mu_y(k)          = mu_y(k) + (j)*glcms(i,j,k);
        end
    end
    mu_x2(k) = mean(p_x(:,k));
    mu_y2(k) = mean(p_y(:,k));
    sigma_x2(k) = std(p_x(:,k));
    sigma_y2(k) = std(p_y(:,k));
end

for k = 1:Nglcms
    for i = 1:L
        for j = 1:L
            sigma_x(k)  = sigma_x(k)  + (((i) - mu_x(k))^2) * glcms(i,j,k);
            sigma_y(k)  = sigma_y(k)  + (((j) - mu_y(k))^2) * glcms(i,j,k);
        end
    end
	sigma_x(k) = sigma_x(k) ^ 0.5; % stddev
	sigma_y(k) = sigma_y(k) ^ 0.5;
end

% Initializing features
f_energ = zeros(Nglcms,1); % Angular Second Moment
f_contr = zeros(Nglcms,1); % Contrast
f_correl = zeros(Nglcms,1); % Correlation
f_sumSqVar = zeros(Nglcms,1); % Sum of squares: Variance
f_invDiffMom = zeros(Nglcms,1); % Inverse Difference Moment
f_sumAvg = zeros(Nglcms,1); % Sum Average
f_sumVar = zeros(Nglcms,1); % Sum Variance
f_sumEntr = zeros(Nglcms,1); % Sum Entropy
f_entr = zeros(Nglcms,1); % Entropy
f_diffVar = zeros(Nglcms,1); % Difference Variance
f_diffEntr = zeros(Nglcms,1); % Difference Entropy
f_infMeasCor1 = zeros(Nglcms,1); % Information Measures of Correlation 1
f_infMeasCor2 = zeros(Nglcms,1); % Information Measures of Correlation 2
if (compMaxCorrCoef)
    f_maxCorCoeff = zeros(Nglcms,1); % Maximal Correlation Coefficient
end
% Additional ones
f_autoCorr = zeros(Nglcms,1); % Autocorrelation
f_dissim = zeros(Nglcms,1); % Dissimilarity
f_clustShade = zeros(Nglcms,1); % Cluster Shade
f_clustProm = zeros(Nglcms,1); % Cluster Prominence
f_maxProb = zeros(Nglcms,1); % Maximum Probability

for k = 1:Nglcms
    for i = 1:L
        for j = 1:L
            
            f_energ(k) = f_energ(k) + (glcms(i,j,k).^2);
			
			f_contr(k) = f_contr(k) + (abs (i - j)) ^ 2 * glcms(i,j,k);
			
            f_entr(k) = f_entr(k) - ...
                (glcms(i,j,k) * log(glcms(i,j,k) + eps));
			
            f_sumSqVar(k) = f_sumSqVar(k) + ((i - glcmMean(k))^2) ...
                * glcms(i,j,k);
			
			f_invDiffMom(k) = f_invDiffMom(k) + (glcms(i, j, k) / ( 1 + (i - j)^2));

			
		end
	end
end

for k = 1:Nglcms
    for i = 1:L
        for j = 1:L
            
			f_autoCorr(k) = f_autoCorr(k) + (i * j * glcms(i,j,k));
			
			f_dissim(k) = f_dissim(k) + (abs(i - j) * glcms(i,j,k));
            
			f_clustShade(k) = f_clustShade(k) + ((i + j - mu_x(k) - mu_y(k)) ^ 3) * glcms(i,j,k);
			
			f_clustProm(k) = f_clustProm(k) + (((i + j - mu_x(k) - mu_y(k)) ^ 4) * glcms(i,j,k));
        end
    end
	
    f_correl(k) = (f_autoCorr(k) - mu_x(k) * mu_y(k)) / (sigma_x(k) * sigma_y(k));
	
	f_maxProb(k) = max(max(glcms(:, :, k)));
end

% Computing Sum Average/Variance/Entropy
for k = 1:(Nglcms)
    for i = 1 : (2 * L - 1)
		
        f_sumAvg(k) = f_sumAvg(k) + (i + 1) * p_xplusy(i, k);
		% f_8 in Haralick1973
        f_sumEntr(k) = f_sumEntr(k) - (p_xplusy(i, k) * log(p_xplusy(i, k) + eps));
    end
	for i = 1 : (2 * L - 1)
		
		f_sumVar(k) = f_sumVar(k) + ((i + 1 - f_sumEntr(k)) ^ 2) * p_xplusy(i,k);
    end
end

% Compute Difference Variance/Entropy 
for k = 1:Nglcms
	for i = 0:(L - 1)
        
		f_diffVar(k) = f_diffVar(k) + (i ^ 2) * p_xminusy(i + 1,k);
		
		f_diffEntr(k) = f_diffEntr(k) - (p_xminusy(i + 1,k) * log(p_xminusy(i + 1,k) + eps));
    end
end

% Compute Information Measures of Correlation and Maximal Correlation
% Coefficient
for k = 1:Nglcms
    hxy(k) = f_entr(k);
    for i = 1:L
        for j = 1:L
            hxy_1(k) = hxy_1(k) - (glcms(i, j, k) * log(p_x(i, k) * p_y(j, k) + eps));
            hxy_2(k) = hxy_2(k) - (p_x(i, k) * p_y(j, k) * log(p_x(i, k) * p_y(j, k) + eps));
			if (compMaxCorrCoef)
                
                for kQ = 1:L
                    Q(i, j, k) = Q(i, j, k) + ...
                      (glcms (i, kQ, k)*glcms(j, kQ, k) / (p_x(i, k)*p_y(kQ, k))); 
                end
            end
        end
        hx(k) = hx(k) - (p_x(i, k) * log(p_x(i, k) + eps));
        hy(k) = hy(k) - (p_y(i, k) * log(p_y(i, k) + eps));
    end

    f_infMeasCor1(k) = (hxy(k) - hxy_1(k)) / (max([hx(k) hy(k)]));

	f_infMeasCor2(k) = (1 - exp(-2 * (hxy_2(k) - hxy(k))))^0.5;
    if (compMaxCorrCoef)
        
        eigenValuesQ(k, :)   = eig(Q (:,:,k)); % Compute Eigen values
        sortedEigenValues(k, :) = sort(eigenValuesQ(k, :), 'descend'); % sort them
        f_maxCorCoeff(k) = sortedEigenValues(k, 2) ^ 0.5;
    end
end

% The order is the same as in Haralick1973 plus the additional aforementioned features
if (compMaxCorrCoef)
    features = [f_energ, f_contr, f_correl, f_sumSqVar, ...
    f_invDiffMom, f_sumAvg, f_sumVar, f_sumEntr, f_entr, ...
    f_diffVar, f_diffEntr, f_infMeasCor1, f_infMeasCor2, ...
    f_maxCorCoeff, f_autoCorr, f_dissim, f_clustShade, ...
    f_clustProm, f_maxProb ...
    ];
else
    features = [f_energ, f_contr, f_correl, f_sumSqVar, ...
    f_invDiffMom, f_sumAvg, f_sumVar, f_sumEntr, f_entr, ...
    f_diffVar, f_diffEntr, f_infMeasCor1, f_infMeasCor2, ...
    f_autoCorr, f_dissim, f_clustShade, ...
    f_clustProm, f_maxProb ...
    ];
end



