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

%% Signal definition
t=0:255;fs=1;
x1 = 2*cos(0.1*pi*t+0.0017*pi*t.^2);
x3 = 2*exp(-0.0015*(t-170).^2).*cos(0.1*pi*t);
x4 = 2*exp(-0.0015*(t-30).^2).*cos(0.5*pi*t);
Signal=x1+x3+x4;



%% QTFD Generation and display
FONT_TYPE='Times-Roman';
FONT_SIZE=14;
LINE_WIDTH=1.5;
tfr_WVD = quadtfd( x1+x3+x4, 255, 1, 'wvd');
t = 0:length(Signal)-1;
f = linspace(0,fs/2,256);
tfr_WVD(tfr_WVD<0)=0;
figure;imagesc(f,t,tfr_WVD.');axis xy; xlabel('Time (sec)');ylabel('Frequency (Hz)');colormap(flipud(hot))
set(gca, 'FontSize',12);
ylabel('Time (s)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
xlabel('Frequency (Hz)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);

%% TSED Generation and Display

[TFD,f]=awvd(Signal,256);
TFD(TFD<0)=0;
t = 0:length(Signal)-1;
nonLinearPlot(TFD.',f.',t);
set(gca, 'FontSize',10);
colormap(flipud(hot));
ylabel('Time (s)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
xlabel('Approximate Frequency (Hz)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
