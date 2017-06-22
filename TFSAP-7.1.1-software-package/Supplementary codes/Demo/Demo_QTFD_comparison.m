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
fs = 32; t = 0:1/fs:8-1/fs;
N = length(t); f = linspace(0,fs/2,256);
s1 = 1*cos(-(6.3578e-09)*(fs^4)*pi*(t-4).^4+0.1*pi*t*fs);
s3 = 1*cos(-(1.734e-08)*(fs^4)*pi*(t-4).^4+0.4*pi*t*fs);
s4 = 1*cos(-(2.7248e-08)*(fs^4)*pi*(t-4).^4+0.6*pi*t*fs);
sig = 1.1*s1 + 0.7*s3 + 0.5*s4;


%% Extended Modified-B Distribution
TFD_EMBD = quadtfd(sig,N-1,1,'emb',0.075,0.5);

%% Spectrogram
[a,f_spec,t_spec,TFD_spec] = spectrogram(sig,hamming(85),84,512,32);
%% Compact kernel Distribution (CKD)
am1 = wvd1(sig,N); A = ckdKernel(1,0.27,0.08,N,N); am = am1.*A;
TFD_CKD = (fft(ifftshift(am,1), [], 1));
TFD_CKD = ifft(fftshift(TFD_CKD,2), [], 2);
TFD_CKD = real(TFD_CKD); TFD_CKD(TFD_CKD<0)=0;

%% S-Method
TFD_SM = specSM(sig,5,85,'hamm',84,512);
%% Plotting
figure; imagesc(f,t,abs(TFD_EMBD')) ;axis xy; set(gca, 'FontSize',12);
xlabel('Frequency (Hz)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
ylabel('Time (s)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
colormap(flipud(hot))

figure; imagesc(f_spec,t_spec,abs(TFD_spec')) ;axis xy;set(gca, 'FontSize',12);
xlabel('Frequency (Hz)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
ylabel('Time (s)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
colormap(flipud(hot))

figure; imagesc(f,t,abs(TFD_CKD'))  ;axis xy;set(gca, 'FontSize',12);
xlabel('Frequency (Hz)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
ylabel('Time (s)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
colormap(flipud(hot))

figure; imagesc(f_spec,t_spec,abs(TFD_SM'))   ;axis xy;set(gca, 'FontSize',12);
xlabel('Frequency (Hz)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
ylabel('Time (s)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
colormap(flipud(hot))