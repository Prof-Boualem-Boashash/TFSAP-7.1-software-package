function featureExtractionProcess(kernel)
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
% This script extract various features from the computed TFDs and organize
% them using three outputs.
% length_features: contains the nulber of epochs for each patient
% signal_class: contains the class of each epoch for each patient; for each
% line one has to take into consideration from epoch 1 to epoch
% length_features(number of the patient) the remaining columns are padding
% with zeros; the number of columns correspond to the maximum number of
% epochs for one epoch(patient 17 in this database).
% signal_feature: this contains the features for each epoch. One line
% represents one patient and as before one has to take into consideration from epoch 1 to epoch
% length_features(number of the patient) the remaining columns are padding
% with zeros;

% clear all; close all
%% PATH
cur_dir = pwd;
% kernel = 'SPEC';
EEGDataSetDir = [cur_dir '\EEGDataSet\']; % DB directory
EEGDataSetFilteredDir  = [EEGDataSetDir 'TFDs\' kernel]; % Filtered  DB directoryocated
addpath(genpath([cur_dir '\code']))
%% PARAMETERS
Np=36; %number of patients
N = 256; %number of samples by epoch
Ne_max = 986; % Maximum number of epochs for one patient
N_features = 89; % Number of feature
signal_class=zeros(Np,Ne_max);
signal_features=zeros(Np,Ne_max,N_features);
length_features=zeros(1,Np);
%% Code for feature extraction on patient by patient basis


for i_patient=1:Np
    
    cd(EEGDataSetFilteredDir)
    eval(['load(''' kernel 'preprocessed_TFD_' num2str(i_patient) '.mat'', ''number_epochs''' ');'])
    
    for i_epoch=1:number_epochs
        eval(['load(''' kernel 'preprocessed_TFD_' num2str(i_patient) '.mat'', ''' kernel '_epoch_' num2str(i_epoch) ''',' '''eegStateN'',''sig_channel_averaged'' );'])
        eval(['tfd=' kernel '_epoch_' num2str(i_epoch) ';']);
        eval(['clear ' kernel '_epoch_' num2str(i_epoch) ';']);
        sig_start = N * (i_epoch-1) + 1;
        sig_end = N * (i_epoch-1) + N;
        sig_cur_epoch = sig_channel_averaged(sig_start:sig_end);
        
        ExtractedFeatures=featureComputation(tfd,sig_cur_epoch);
        text = sprintf('patient: %d, segment: %d .',i_patient,i_epoch);
        disp(text);
        
        signal_class(i_patient,i_epoch)=eegStateN(i_epoch,3);
        signal_features(i_patient,i_epoch,1:length(ExtractedFeatures))=ExtractedFeatures;
    end
    length_features(i_patient)=number_epochs;
end
cd(cur_dir);
% path = [cur_dir '\' kernel 'results_feature_extraction'];
eval(['save(''result_' kernel ''',''signal_class'',''signal_features'',''length_features'')']);