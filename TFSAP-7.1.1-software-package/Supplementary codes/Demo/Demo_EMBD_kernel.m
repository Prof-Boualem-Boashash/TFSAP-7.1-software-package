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
% Description:
%   This funcion allows the reader to reproduce Fig. 6 of the manuscript
%   representing two different EMBD kernel in the ambiguity domain.

%% Parameters
alpha_1 =0.1; beta_1 = 0.2;win_N = 256;
alpha_2 =0.3; beta_2 = 0.1;
FONT_TYPE='Times-Roman';
FONT_SIZE=14;
LINE_WIDTH=1.5;

g_embd_1 = extnd_mbd(alpha_1,beta_1,0.5,win_N);
g_embd_2 = extnd_mbd(alpha_2,beta_2,0.5,win_N);

tau = linspace(-win_N/2,win_N/2,win_N);
Doppler = linspace(-0.5,0.5,win_N);

figure;imagesc(tau,Doppler,g_embd_1);axis xy;set(gca, 'FontSize',12);
xlabel('Lag (s)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
ylabel('Doppler (Hz)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
colormap(flipud(hot)); 

figure;imagesc(tau,Doppler,g_embd_2);axis xy;set(gca, 'FontSize',12);
xlabel('Lag (s)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
ylabel('Doppler (Hz)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
colormap(flipud(hot));