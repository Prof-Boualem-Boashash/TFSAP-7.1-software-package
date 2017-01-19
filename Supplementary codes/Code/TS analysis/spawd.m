
function TFD = spawd(x, nh0, ng0, M)
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
% stationary signals: a tutorial review with features performance comparison
%, Digital Signal Processing, In Press.
% [2] B. Boashash, Samir Ouelha, Efficient Software Matlab Package to compute
% Time-Frequency Distributions and related Time-Scale methods with extraction
% of signal characteristics, SoftwareX, In Press.
%
% This study was funded by grants from the ARC and QNRF NPRP 6-885-2-364
% and NPRP 4-1303-2-517 
%
% Last Modification: 25-12-2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
%
% Computes the smoothed_pseudo Affine Wigner-Ville Distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Computed Variables
thr  = 1e-2;
fmin = thr;
fmax = 0.5;
N = length(x);
q = (fmax/fmin)^(1/(M-1));
a = exp(((1:M)-1).*log(q));
f_geo = fmin*a;

%% Time Smoothing Window G
umax  = log(fmax/fmin);
Teq   = nh0/(fmax*umax);
if(Teq < 2*nh0), M0 = (2*nh0^2)/Teq-nh0+1; else M0 = 0; end;
MU    = round(nh0+M0);
if(ng0 == 0)
    G = ones(2*MU,1);
else
    sigma_t = ng0*fmax/sqrt(6*log(10));
    a_u = 2 * pi^2 * sigma_t^2 * umax^2 / log(10) ;
    G = (exp(-(a_u*log(10)/MU^2)*(-MU:MU-1).^2))';
end

%%
umin = -umax;
u    = linspace(umin,umax,2*MU+1);
u    = u(1:2*MU);
u(MU+1) = 0;
beta = ((0:(2*M-1))/M-1)./(2*log(q));
l1 = zeros(2*MU,2*M);
l2 = zeros(2*MU,2*M);
for m = 1:2*MU
    temp = u(m);
    y = ones(size(u(m)));
    ind = find(temp ~= 0);
    y(ind) = 2*(exp(-temp(ind))-1)./(exp(-2*temp(ind))-1);
    l1(m,:) = exp(-1i*2*pi*beta*log(y));
    y(ind) = 2*(exp(temp(ind))-1)./(exp(2*temp(ind))-1);
    l2(m,:) = exp(-1i*2*pi*beta*log(y));
end

%% Wavelet Decomposition
s = real(x);
matxte = zeros(M,N);
[~,wt] = scalogram(s,nh0, M);
for ptr = 1:M
    matxte(ptr,:) = wt(ptr,:).*sqrt(a(M-ptr+1));
end

%% Computing the SPAWD
TFD = zeros(M,N);
S  = zeros(1,2*M);
for n = 1:N
    S(1:M) = matxte(:,n).';
    Mellin = fftshift(ifft(S));
    MX1 = fft((l1.*repmat(Mellin,2*MU,1)).');
    MX2 = fft((l2.*repmat(Mellin,2*MU,1)).'); 
    TX1 = MX1(1:M,:).';
    TX2 = MX2(1:M,:).';    
    waf = real(TX1.*conj(TX2)).*repmat(G,1,M);
    TFD(:,n) = (sum(waf).*f_geo).';
end

end