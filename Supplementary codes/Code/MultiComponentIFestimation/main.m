
clear; close all; clc;

%% Loading the Bat Signal
sig = load ('bat1'); 

%% Parameters
N  = length(sig);        % Number of Samples
M  = 512;                % Frequency bins
fs = 142;                % Sampling Frequency (in KHz)
t  = 0:1/fs:(N-1)/fs;    % Time Array
f  = linspace(0,fs/2,M); % Frequency Array

%% Computing Time-Frequency Distributions
TFD = cmpt(sig, 'ckd', 4, 0.1, 0.1, M);

%% Component Linking Algorithm
thresh    = 0.01;
minlength = 32;
[edgeim, peaks] = component_linking(TFD, thresh, minlength);

%% BSS-STRE Algorithm
np = 4; W = 25; thr1 = 0.1; thr2 = 0.03;
[IF_est, lim, timesupport] = BSS_STRE(TFD, np, W, thr1, thr2);

%% Plotting
tfr = 2;
figure; tfsapl(sig,TFD(:,1:tfr:end),'TimePlot','off','GrayScale','on','SampleFreq',fs,...
    'xlabel','Frequency (kHz)','title','Original Signal TFD');

figure; imagesc(f,t,1-edgeim'); ylabel('Time (s)'); xlabel('Frequency (kHz)');
colormap gray; title('Estimated IF using Component Linking Algorithm'); axis xy;

figure;
for i = 1:np
    plot(fs*IF_est(lim(1,i):lim(2,i),i),t(lim(1,i):lim(2,i)),'k'); hold on;
end
axis([0 f(end) 0 t(end)]); xlabel('Frequency (kHz)');
colormap gray; title('Estimated IF using BSS-STRE Algorithm'); axis xy;
