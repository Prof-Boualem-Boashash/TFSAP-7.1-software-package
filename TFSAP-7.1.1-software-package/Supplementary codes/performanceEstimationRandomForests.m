function [predictedOutput, RealLabel, Performance, ACC, Sen, Spe] = performanceEstimationRandomForests(length_features, signal_features, signal_class, index)
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
% This function estimates the classification performance using random
% forest classifier with 500 trees. The method used for the estimation is
% the leave-one-out-patient.
% predictedOutput represents the estimated output for each epoch.
% RealLabel represents the real outputs for each epoch.
% Performance is the performance for each patient.
% Accuracy is the total classification accuracy.
% Sen is the sensitivity and Spe is the specificity.


predictedOutput = [];
RealLabel = [];
Np = length(length_features);
Performance = zeros(1,Np);
for tt=1:Np
    mask=[];
    TEST=[];
    normal_feature_vector=[];
    artifact_feature_vector=[];
    seizure_feature_vector=[];
    k1=0;
    k2=0;
%     k3=0;
    
    %Data for learning
    for ii=1:Np
        if ii~=tt % one patient out
            
            for jj=1:length_features(ii)
                if  signal_class(ii,jj)==3
                    k1=k1+1;
                    z(1:length(index))=signal_features(ii,jj,index);
                    normal_feature_vector(k1,:)=z;
                elseif signal_class(ii,jj)==4
                    k1=k1+1;
                    z(1:length(index))=signal_features(ii,jj,index);
                    normal_feature_vector(k1,:)=z;
                else
                    k2=k2+1;
                    z(1:length(index))=signal_features(ii,jj,index);
                    seizure_feature_vector(k2,:)=z;
                end
            end
            
        end
    end
    
    Training=[seizure_feature_vector( :,:);normal_feature_vector(:,:);artifact_feature_vector(:,:)];
    Group     = [ones(k2,1); zeros(k1,1)];
    
    for jj=1:(length_features(tt))
        z(1:length(index))=signal_features(tt,jj,index);
        TEST(jj,:)=z;
        if signal_class(tt,jj)==3 %non seizure and artefact
            mask(jj) = 0;
        elseif signal_class(tt,jj)==4
            mask(jj) = 0;
        else
            mask(jj) =  1;
        end
    end
    
    
    %% Training
    rng(1);
    B = TreeBagger(100,Training,Group,'Method', 'classification');
    
    %% Estimation Performance
    predChar1 = B.predict(TEST);
    % Predictions is a cell though. We want it to be a number.
    predictedClass = str2double(predChar1); 
    predictedOutput = [predictedOutput;predictedClass];
    RealLabel = [RealLabel mask];
    Performance(tt)=sum(predictedClass(:)==mask(:))/length(mask);
end

%% Results
ACC = sum(Performance(:).*length_features(:))/sum(length_features);

TP=sum(predictedOutput(:).*RealLabel(:));% Y=1,output=1
FP=sum(predictedOutput(:).*(1-RealLabel(:)));% Y=0,output=1
TN=sum((1-predictedOutput(:)).*(1-RealLabel(:)));%Y=0,output=0
FN=sum((1-predictedOutput(:)).*(RealLabel(:)));%Y=1,output=0

% PPV = TP/(TP+FP);
% NPP = TN/(TN+FN);
Sen=TP/(TP+FN);
Spe=TN/(TN+FP);


% Ref:
% B. Boashash, S. Ouelha, Automatic signal abnormality detection using
% time-frequency features and machine learning: a newborn EEG seizure
% case study., Knowledge-Based Systems. doi:10.1016/j.knosys.2016.05.027.

