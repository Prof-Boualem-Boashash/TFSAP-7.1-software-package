%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% Authors: Prof. Boualem Boashash   (boualem.boashash@gmail.com)
%          Samir Ouelha             (samir_ouelha@hotmail.fr)
%          
% The following references should be cited whenever this script is used:
% [1] B. Boashash, Samir Ouelha, Designing high-resolution time-frequency
% and time-scale distributions for the analysis and  classification of non 
% stationary signals: a tutorial review with a comparison of features performance
%, Digital Signal Processing, In Press, 2017.
% [2] B. Boashash, Samir Ouelha, Efficient Software Matlab Package to compute
% Time-Frequency Distributions and related Time-Scale methods with extraction
% of signal characteristics, SoftwareX, In Press, 2017.
%
% This study was funded by grants from the ARC and QNRF NPRP 6-885-2-364
% and NPRP 4-1303-2-517 
%
%
addpath(genpath('..\Code'));
[Signal,fs] = audioread('sp1.wav');
t = 0:1/fs:(length(Signal)-1)/fs;

%% TSED Generation and Display
FONT_TYPE='Times-Roman';
FONT_SIZE=14;
LINE_WIDTH=1.5;

[a,f_spec,t_spec,TFD_spec] = spectrogram(Signal,hamming(512),500,512,fs);
figure; imagesc(f_spec,t_spec,10*log10(TFD_spec')) ;axis xy;set(gca, 'FontSize',12);
xlabel('Frequency (Hz)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
ylabel('Time (s)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
colormap(jet)

[TFD,CWT,f] = scalogram(Signal,50, 256);
nonLinearPlot(10*log10(TFD'),fs*f',t);
set(gca, 'FontSize',12);
ylabel('Time (s)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
xlabel('Approximate Frequency (Hz)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);colormap(jet)