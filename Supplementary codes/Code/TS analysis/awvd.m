function  TFD = awvd(x, M)

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
%
% Description:
%

% This function generates the Affine Wigner-Ville distribution.

%	x : signal (in time) to be analyzed. 

%	TFD : time-frequency matrix containing the coefficients of the
%	   distribution (x-coordinate corresponds to uniformly sampled
%	   time, and y-coordinate corresponds to a geometrically sampled
%	   frequency). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Computed Variables
M = 2*M;
fmin = 1e-2;
fmax = 0.5;
N = length(x);
MM = (N+rem(N,2))/2;
ff = @(u)(exp(u)-fmax/fmin);
umax = fzero(ff,[0 4]);
Teq = MM/(fmax*umax);
if(Teq < N)
    M0 = round((2*MM^2)/Teq-MM)+1;
    M1 = MM + M0;
elseif(Teq >= N)
    M1 = MM;
end

%% Geometric sampling of the analyzed spectrum
x = x(:);
s = real(x);
q = (fmax/fmin)^(1/(M-1));
t = (1:N)-MM-1;
geo_f = fmin*(exp(((1:M)-1).*log(q)));
tfmatx = exp(-2*1i*t'*geo_f*pi);
S = s'*tfmatx;
S(M+1:M*2) = 0;

%% Mellin transform computation of the analyzed signal
Mellin = fftshift(ifft(S));
du = abs(umax)/(M1);
u(1:2*M1) = -umax:du:umax-du;
u(M1+1) = 0;
beta = ((0:(2*M-1))/M-1)/(2*log(q));

%% Computation of P0(t.f,f)
waf = zeros(2*M1,M);
for n = [1:M1,M1+2:2*M1],
    MX1 = exp((-2*1i*pi*beta+0.5)*log((u(n)/2)*exp(-u(n)/2)/sinh(u(n)/2))).*Mellin;
    MX2 = exp((-2*1i*pi*beta+0.5)*log((u(n)/2)*exp( u(n)/2)/sinh(u(n)/2))).*Mellin;
    FX1 = fft(fftshift(MX1)) ;
    FX1 = FX1(1:M) ;
    FX2 = fft(fftshift(MX2)) ;
    FX2 = FX2(1:M) ;
    waf(n,:) = FX1.*conj(FX2);
end
waf(M1+1,:) = S(1:M).*conj(S(1:M));
waf  = [waf(M1+1:2*M1,:) ; waf(1:M1,:)].*geo_f(ones(2*M1,1),:);
tffr = ifft(waf);
tffr = real(rot90([tffr(M1+1:2*M1,:) ; tffr(1:M1,:)],-1));

%% Conversion from [t.f,f] to [t,f] using a 1-D interpolation
temp  = zeros(M,N);
Ts2   = (N-1)/2 ;
gamma = linspace(-geo_f(M)*Ts2,geo_f(M)*Ts2,2*M1) ;
for n = 1:M
    ind = find(gamma>=-geo_f(n)*Ts2 & gamma<=geo_f(n)*Ts2);
    xx  = gamma(ind);
    yy  = tffr(n,ind);
    xi = ((1:N)-Ts2-1)*geo_f(n);
    v  = interp1(xx,yy,xi,'pchip');
    temp(n,:) = v(:).';
end
TFD = zeros(M/2,N);
for i = 1:N
    TFD(:,i) = downsample(temp(:,i),2);
end

end