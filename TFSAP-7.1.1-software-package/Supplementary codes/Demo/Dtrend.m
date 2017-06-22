function [z_stat,z_trend] = Dtrend(z,lambda)
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

[P,Q]=size(z);
if (P<Q)
    z=z.';
end


T=length(z);
Id=speye(T);
D2=spdiags(ones(T-2,1)*[1 -2 1],0:2, T-2,T);
z_stat=(Id-inv(Id+lambda^2*(D2'*D2)))*z;
z_trend=z-z_stat;

end

