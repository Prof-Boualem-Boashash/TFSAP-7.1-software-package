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
FONT_TYPE='Times-Roman';
FONT_SIZE=14;
LINE_WIDTH=1.5;

%% Signal definition
fs = 32;
n = 0:1/fs:8-1/fs; N = length(n);
sig1 = cos(2*pi*n*1);
sig2 = cos(2*pi*n*3);
sigm = 0.0020;
B = fir1(16*2,0.08,'high');
s1 = 10*exp(-(fs^2*(n-0.9375).^2)/sigm);
s1 = filter(B,1,s1);
s2 = 10*exp(-(fs^2*(n-2.1875).^2)/sigm);%.*cos(0.1*pi*n);
s2 = filter(B,1,s2);
s3 = 10*exp(-(fs^2*(n- 3.4375).^2)/sigm);%.*cos(0.1*pi*n);
s3 = filter(B,1,s3);
s4 = 10*exp(-(fs^2*(n-4.6875).^2)/sigm);%.*cos(0.1*pi*n);
s4 = filter(B,1,s4);
s5 = 10*exp(-(fs^2*(n- 5.9375).^2)/sigm);%.*cos(0.1*pi*n);
s5 = filter(B,1,s5);
s = s1 + s2 + s3 + s4 + s5;
Signal = sig1 + sig2 + s;



%% TFD Generation and display
disp('1: Spectrogram with short window');
[~,f,t,P]=spectrogram(Signal,24,23,256,fs);
figure;imagesc(f,t,P');axis xy;set(gca, 'FontSize',12);
ylabel('Time (s)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
xlabel('Frequency (Hz)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
colormap(flipud(hot))
% title('Spectrogram with short window');

disp('2: Spectrogram with long window');
[~,f,t,P]=spectrogram(Signal,84,83,256,fs);
figure;imagesc(f,t,P');axis xy;set(gca, 'FontSize',12);
ylabel('Time (s)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
xlabel('Frequency (Hz)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
% title('Spectrogram with long window');
colormap(flipud(hot))

