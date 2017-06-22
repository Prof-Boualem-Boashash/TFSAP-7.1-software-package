% 
% *Copyright 2016 Boualem Boashash
% 
%  Licensed under the Apache License, Version 2.0 (the "License");
%  you may not use this file except in compliance with the License.
%  You may obtain a copy of the License at
% 
%      http://www.apache.org/licenses/LICENSE-2.0
% 
%  Unless required by applicable law or agreed to in writing, software
%  distributed under the License is distributed on an "AS IS" BASIS,
%  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%  See the License for the specific language governing permissions and
%  limitations under the License.
% 
%  Author:                 Boualem Boashash         (boualem.boashash@gmail.com)
%  Maintainer since 2015:  Samir Ouelha  			(samir_ouelha@hotmail.fr)
% 
%  The following references 2 should be cited whenever this script is used:
%  [1] B. Boashash, Samir Ouelha, Designing time-frequency  and time-scale
%  features for efficient classification of non stationary signals: a tutorial
%  review with performance comparison, Digital Signal Processing, In Press.
%  [2] B. Boashash, Samir Ouelha, Efficient Software Matlab Package for the calculation
%  of Time-Frequency Distributions related Time-Scale methods and the extraction of
%  signal characteristics, SoftwareX, In Press.
%  In addition, the following 3rd reference is relevant for use of the executable
%  [3] B. Boashash (ed.), Time-Frequency Signal Analysis and Processing, 2nd Edition
%  (London: Elsevier / Academic Press, December 2015); ISBN 978-0-12-398499-9. Book website:
%  http://booksite.elsevier.com/9780123984999/
%  Lastly the following reference is useful for understanding the basics of Instantaneous
%  Frequency estimation:
%  [4] B. Boashash, "Estimating and interpreting the instantaneous frequency of
%  a signal—part 2: algorithms and applications", Proc. IEEE 80 (4) (1992) 540-568.
% 
%  Description:
% 
% Multi-Directional Distribution

function [tfdF,A,am1] = mdd(sig,c,D,E,theta,tr)

D = D*0.1;
sig=sig(:)';
am1 =wvd1(sig);

%% processing
A=kernel_MDD(sig,c,D,E, theta);
am=am1.*A;
tfd = (fft(ifftshift(am,1), [], 1));
tfd=  ifft(fftshift(tfd,2), [], 2);
tfd=real(tfd);
tfd(tfd<0)=0;
tfdF = tfd(:,1:tr:end);
end

function [amb, tfrep] = wvd1(x)

N = length(x);
analytic_sig_ker = signal_kernal(x);
tfrep = real(1./N.*fft(ifftshift(analytic_sig_ker,1), N, 1));
amb = fftshift(1./N.*fft(analytic_sig_ker, N, 2),2);
end


function analytic_sig_ker = signal_kernal(x)
N = length(x);
if isreal(x)
    if mod(length(x),2) == 0
        true_X = fft(x);
        analytic_X = [true_X(1) 2.*true_X(2:N/2) true_X(N/2+1) zeros(1,N/2-1)];
        analytic_x = ifft(analytic_X);
    else
        true_X = fft(x);
        analytic_X = [true_X(1) 2.*true_X(2:ceil(N/2)) zeros(1,floor(N/2))];
        analytic_x = ifft(analytic_X);
    end
else
    analytic_x=x;
end
analytic_sig_ker = zeros(N,N);

for m = -round(N/2-1):1:round(N/2-1)
    analytic_sig_ker(m+round(N/2)+1,:) = sig_ker_corr(analytic_x,m);
end
end


function sig_ker_m = sig_ker_corr(x,m)

N = length(x);
z_nmm = [zeros(1,N+m) x zeros(1,N-m)];
z_npm = [zeros(1,N-m) x zeros(1,N+m)];
ker_nm = z_npm.*conj(z_nmm);
sig_ker_m = ker_nm(N+1:2*N);
end

function A=kernel_MDD(sig,c,D,E, theta)
%%%% sig : represent signal
%%% theta : vector of angles
%%% c : paramter to adjust the form of shape
%%% D : paramter to adjust the size of window
%%% E : paramter to adjust the size of window

theta=(theta)*pi/180;
N1=length(sig);
if mod(N1,2)==0;
    N=N1;
else
    N=N1;
    N1=N1+1;
end
y=(-.5+1/N1:1/N1:.5-1/N1);
x=(-.5+1/N:1/N:.5-1/N);
A=zeros(N1,N);
for kk=1:length(theta)
    for jj=1:N1-1
        for ii=1:N-1
            x1=cos(theta(kk))*x(ii)-sin(theta(kk))*y(jj);
            y1=sin(theta(kk))*x(ii)+cos(theta(kk))*y(jj);
            
            z1=(exp(c(kk)*(1-exp(abs(x1/D(kk)).^2))));
            w=(exp(c(kk)*(1-exp(abs(y1/E(kk)).^2))));
            z1=z1*w;
            if and(abs(x1)<D(kk),abs(y1)<E(kk))
                A(ii,jj)=A(ii,jj)+z1;
            end
        end
    end
end

A=A/norm(A);
end
