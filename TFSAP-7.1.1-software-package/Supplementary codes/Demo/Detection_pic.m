function [ peak,loc ] = Detection_pic(sig)

%  Author:
%     Samir Ouelha, Postdoc for Prof. Boashash
%
%  Last update:
%     21/04/2016
%Please cite the following reference:

%   [1] B.Boashash and S. Ouelha
%       “Improved design of high-resolution Quadratic Time-Frequency
%       distributions for the analysis of multicomponent non-stationary signals”
%       IEEE Transactions on Signal Processing, 2016
%
% This work was supported by Qatar Foundation grant NPRP 4-1303-2-517.

sig = sig(:)';
sa=sign(diff([-inf sig]));

sb=sign(diff([-inf sig(end:-1:1)]));
sb=sb(end:-1:1);

idx=(sa==1 & sb==1);

loc=find(idx==1);
peak=sig(loc);

end