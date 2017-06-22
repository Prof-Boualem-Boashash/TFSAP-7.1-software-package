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
%  The following 2 references should be cited whenever this script is used:
%  [1] B. Boashash, Samir Ouelha, Designing time-frequency  and time-scale
%  features for efficient classification of non stationary signals: a tutorial
%  review with performance comparison, Digital Signal Processing, 2017.
%  [2] B. Boashash, Samir Ouelha, Efficient Software Matlab Package for the calculation
%  of Time-Frequency Distributions related Time-Scale methods and the extraction of
%  signal characteristics, SoftwareX, 2017.
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
%  Wigner-Distribution
function tfdF = wd(x,fftl,tr)


sig_ker = signal_kernal(x);
tfd = real(fft(ifftshift(sig_ker,1), fftl, 1));
tfdF = tfd(:,1:tr:end);

end

function sig_ker = signal_kernal(x)
N = length(x);
sig_ker = zeros(N,N);
for m = -round(N/2-1):1:round(N/2-1)
    sig_ker(m+round(N/2)+1,:) = sig_ker_corr(x,m);
end
end


function sig_ker_m = sig_ker_corr(x,m)

N = length(x);
z_nmm = [zeros(1,N+m) x zeros(1,N-m)];
z_npm = [zeros(1,N-m) x zeros(1,N+m)];
ker_nm = z_npm.*conj(z_nmm);
sig_ker_m = ker_nm(N+1:2*N);
end

