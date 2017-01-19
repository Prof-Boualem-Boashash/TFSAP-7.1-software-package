function [ACC, Sen, Spe,C,Performance] = main(kernel, featureSelection)


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
% Description:
% This is the main function for generating the classification results. We
% can choose the TFD kernel, the method for filter feature selection.
% After this, one can use the feature selection using the wrapper methods
% combined with SFFS method (it takes long time).
%
%

%% Parameters

nbFeatures = 89; %total number of features
nbFeatureFiltered = 80; %total number conserved after filter feature selection
%% PATH
cur_dir = pwd;
EEGDataSetDir = [cur_dir '\EEGDataSet\']; % DB directory
pathTFD = [EEGDataSetDir 'TFDs\' kernel];
pathPreprocess = [EEGDataSetDir 'Filtered'];
addpath(genpath([cur_dir '\code']));
%replace "path_to_TFSAP_ToolBox" with the actual path to the toolbox on
%your PC
%addpath('path_to_TFSAP_ToolBox');

%% Pre processing of the original data;
ext = '*.mat';
files = dir(fullfile(pathPreprocess,ext));
if  exist(pathPreprocess, 'file')==7
    if isempty(files)
        disp([pathPreprocess,' is empty']);
        PreProcessEEGdata;
    else
        disp('The signals are computed');
    end
else
    mkdir(pathPreprocess);
    PreProcessEEGdata;
end
%% Computation of the TFDs

ext = '*.mat';
files2 = dir(fullfile(pathTFD,ext));
if  exist(pathTFD, 'file')==7
    if isempty(files2)
        disp([pathTFD,' is empty']);
         computeTFDs(kernel);
    else
        disp('The TFDs are computed; there is no need to recompute it again');
    end
else
    mkdir(pathTFD);
    computeTFDs(kernel);
end

%% Feature extraction
if  exist([cur_dir '\result_' kernel '.mat'], 'file')~=2
    featureExtractionProcess(kernel);
end
load([cur_dir '\result_' kernel]);

%% Feature selection process
if featureSelection ==true
    disp('Starting of the feature selection process');
    [Data, Label]=dataManager(length_features, signal_features, signal_class, 1:nbFeatures);
    for k=1:length(Label)
        if Label(k)==-1 %non seizure and artefact
            Label(k) = 1;
        elseif Label(k)==0
            Label(k) = 1;
        else
            Label(k) =  2;
        end
    end
    [~,IX]  = maximalMarginalDiversity(Data(:,1:nbFeatures),Label);
    signal_features_reduced = signal_features(:,:,IX(1:nbFeatureFiltered));
    save(['result_' kernel],'signal_features_reduced','-append');
    [Xbest, J_hat_best]=SFFS_test(nbFeatureFiltered,kernel);
    disp('End of the feature selection process');
end
%% Performance Estimation using random Forest
rng(1);
[predictedOutput, RealLabel, Performance, ACC, Sen, Spe] = performanceEstimationRandomForests(length_features, signal_features, signal_class, 1:nbFeatures);
[C, order] = confusionmat(RealLabel(:),predictedOutput(:));disp(order);
disp(['The accuracy of the ' kernel ' kernel is: ' num2str(100*ACC) '%'])
disp(['The sensitivity of the ' kernel ' ketnel is: ' num2str(100*Sen) '%'])
disp(['The specificity of the ' kernel ' ketnel is: ' num2str(100*Spe) '%'])

% Ref:
% B. Boashash, S. Ouelha, Automatic signal abnormality detection using
% time-frequency features and machine learning: a newborn EEG seizure
% case study., Knowledge-Based Systems. doi:10.1016/j.knosys.2016.05.027.



