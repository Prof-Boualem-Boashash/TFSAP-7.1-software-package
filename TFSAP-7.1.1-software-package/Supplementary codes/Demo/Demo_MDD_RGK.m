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
clear all; close all;
FONT_TYPE='Times-Roman';
FONT_SIZE=14;
LINE_WIDTH=1.5;
%% Signal definition and parameters
t=0:255;  % time vector
f = linspace(0,0.5,length(t)); % frequency vector
fs=1;     % Sampling frequency
th=5;     % Threshold for the display
tr=5;     % time resolution
nbSig = 4;

s1 = 1* cos(2*pi*(0.0006)*t.*t+2*pi*0.15*t).*[zeros(1,20) ones(1,180) zeros(1,56)];
s2 = sin(2*pi*0.03*t).*[zeros(1,20) ones(1,180) zeros(1,56)];
s3 = sin(2*pi*0.15*t).*[zeros(1,80) ones(1,150) zeros(1,26)];

sig = s1 + s2  + s3;

%% Multi-directional distribution
[TFD_MDD,A,am1] = myMDD(sig(:)',1,0.2,0,0);

%% Radially Gaussian Kernel
[TFD_RGK, Phi] = rgk(sig,2);

%% Plotting
figure; imagesc(f,t,abs(TFD_MDD')) ;axis xy; set(gca, 'FontSize',12);
xlabel('Frequency (Hz)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
ylabel('Time (s)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
colormap(flipud(hot))

figure; imagesc(f,t,abs(TFD_RGK')) ;axis xy;set(gca, 'FontSize',12);
xlabel('Frequency (Hz)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
ylabel('Time (s)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
colormap(flipud(hot))