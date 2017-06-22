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
%%
load EEG_examples
fs = 32; t = 0:1/fs:8-1/fs;
N = length(t); f = linspace(0,fs/2,N);
FONT_TYPE='Times-Roman';
FONT_SIZE=14;
LINE_WIDTH=1.5;

TFD_EMBD_seiz = quadtfd(sig_seiz,N-1,1,'emb',0.075,0.5);
LBP_seiz = myLBP(TFD_EMBD_seiz);

TFD_EMBD_back = quadtfd(sig_back,N-1,1,'emb',0.075,0.5);
LBP_back = myLBP(TFD_EMBD_back);



figure; imagesc(f,t,abs(TFD_EMBD_seiz')) ;axis xy; set(gca, 'FontSize',12);
xlabel('Frequency (Hz)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
ylabel('Time (s)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
colormap(flipud(hot))

figure; imagesc(f,t,abs(TFD_EMBD_back')) ;axis xy; set(gca, 'FontSize',12);
xlabel('Frequency (Hz)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
ylabel('Time (s)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
colormap(flipud(hot))

figure; imagesc(f,t,abs(LBP_seiz')) ;axis xy; set(gca, 'FontSize',12);
xlabel('Frequency (Hz)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
ylabel('Time (s)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
colormap(jet)

figure; imagesc(f,t,abs(LBP_back')) ;axis xy; set(gca, 'FontSize',12);
xlabel('Frequency (Hz)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
ylabel('Time (s)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
colormap(jet)