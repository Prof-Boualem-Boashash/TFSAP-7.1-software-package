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
%

function [g_extmb] = extnd_mbd(a,b,min_fre_diff,win_N)

if nargin == 3
    win_N = 128;
end
if nargin == 2
    win_N = 128; min_fre_diff = 0.5;
end

G_mb_Dopper = zeros(win_N,win_N); 
for n = -win_N/2:win_N/2
    G_mb_Dopper(:,mod(n,win_N)+1) = cosh( n ).^( -2 * a );
end
tmp1 = fft(G_mb_Dopper(1,:));
G_mb_Dopper = G_mb_Dopper./tmp1(1);
g_mb_Dopper = real(fft(G_mb_Dopper.').');

G_mb_lag = zeros(win_N,win_N); 
for n = -win_N/2:win_N/2
    G_mb_lag(:,mod(n,win_N)+1) = cosh( n ).^( -2 * b );
end
tmp1 = fft(G_mb_lag(1,:));
G_mb_lag = G_mb_lag./tmp1(1);
g_mb_lag = real(fft(G_mb_lag.').');

effective_win = floor(min_fre_diff*win_N);
g_extmb = g_mb_lag.'.*g_mb_Dopper;
tmp = -win_N/2:win_N/2-1;
g_extmb(:,abs(tmp) <= (win_N/2-effective_win)) = 0;

g_extmb = fftshift(g_extmb);


