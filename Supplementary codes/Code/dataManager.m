function [X, Group]=dataManager(length_features_spec, signal_features_spec, signal_class_spec, index)
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
% This function organize the data in another way to be able to use filter
% feature selection algorithms.
% X is composed of all features for all patients; the last column is the
% number of the patient.
% Group is the label for each epoch.


Np = length(length_features_spec);
normal_feature_vector=[];
artifact_feature_vector=[];
seizure_feature_vector=[];
k1=0;
k2=0;
k3=0;

%Data for learning
for ii=1:Np
    
    for jj=1:length_features_spec(ii)
        
        if  signal_class_spec(ii,jj)==3
            k1=k1+1;
            z(1:length(index))=signal_features_spec(ii,jj,index);
            normal_feature_vector(k1,:)=[z ii];
        elseif signal_class_spec(ii,jj)==4
            k2=k2+1;
            z(1:length(index))=signal_features_spec(ii,jj,index);
            artifact_feature_vector(k2,:)=[z ii];
        else
            k3=k3+1;
            z(1:length(index))=signal_features_spec(ii,jj,index);
            seizure_feature_vector(k3,:)=[z ii];
        end
        
    end 
    
end


X =[seizure_feature_vector( :,:);normal_feature_vector(:,:);artifact_feature_vector(:,:)];
Group     = [ones(k3,1); zeros(k1,1);-ones(k2,1)];

