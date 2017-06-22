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
% This function preprocesses the original signal filtering them and
% subsampling them.


% PreProcessBrisbEEGDB


%% PATH
cur_dir = pwd;

%Code for IIR filtering
EEGDataSetDir = ['EEGDataSet\']; % DB directory
% EEGDataSetArtRemSubDir = ['artifact_removed\multi\'];
% EEGDataSetDirArtRem = [EEGDataSetDir EEGDataSetArtRemSubDir]; %

%% PARAMETERS
fs = 256; %sampling frequency
lfc = 0.8; hfc = 12; % bandpass filter cut-off frequencies
T = 8; % epoch time in seconds
N = T*fs; % epoch samples
channels = 1:20;
Np=36; %number of patients

for ii=1:36
    filenameN=['original\' num2str(ii) '_data.mat'];
    pathN = [EEGDataSetDir filenameN];
    XN = load(pathN);
    dataN = XN.data; % eeg data
    maskN = XN.mask; % seizure mask per neurologist (in seconds)
    eegStateN = XN.eegState; % contains seizure event number (0: for background), epoch # and class
    fdataN = iirfilt(dataN',fs,lfc,0).';
    fdataN = iirfilt(fdataN.',fs,0,hfc).'; % bandpass filtering
    sig_multichannel=resample(fdataN,32,256);
    sig_channel_averaged = mean(sig_multichannel,2);
    pathNew = [EEGDataSetDir 'Filtered\preprocessed_data_' num2str(ii)];
    save(pathNew,'sig_multichannel', 'dataN', 'eegStateN', 'sig_channel_averaged'); 
    % This file contains the preprocessed multi-channel signals (fs = 32
    % Hz), the averaged corresponding signal, the original non processed
    % signals, and the original EEG state for further classification
    % (label).
end

