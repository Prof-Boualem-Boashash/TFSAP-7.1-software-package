function A=kernel_MDD(sig,c,D,E, theta)
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


%%%% sig : represent signal
%%% theta : vector of angles
%%% c : paramter to adjust the form of shape
%%% D : paramter to adjust the size of window
%%% E : paramter to adjust the size of window



theta=(theta)*pi/180;


N1=length(sig);N=N1;
y=(-.5+1/N1:1/N1:.5-1/N1);
x=(-.5+1/N:1/N:.5-1/N);
A=zeros(N,N);
for kk=1:length(theta)
    for jj=1:N1-1
        for ii=1:N-1
            x1=cos(theta(kk))*x(ii)-sin(theta(kk))*y(jj);
            y1=sin(theta(kk))*x(ii)+cos(theta(kk))*y(jj);
            
            z1=exp(c*(1-exp(abs(x1/D(kk)).^2)));
            w=exp(c*(1-exp(abs(y1/E(kk)).^2)));
            z1=z1*w;
            if and(abs(x1)<D(kk),abs(y1)<E(kk))
                A(ii,jj)=A(ii,jj)+z1;
            end
        end
    end
end

A=A/max(A(:));
% for kk=1:length(theta)
%     for ii=1:length(y)
%         for jj=1:length(x)
%             x1=cos(theta(kk))*x-sin(theta(kk))*y;
%             y1=sin(theta(kk))*x+cos(theta(kk))*y;
%             chk_val = (x1(jj)^2)/(D(kk)^2)+ (y1(jj)^2)/(E(kk)^2);
%             if (abs(x1(jj))<abs(D(kk)) && abs(y1(ii))<abs(E(kk)))
% %             if(chk_val<1)%(x(jj).^2/E.^2+ y(ii).^2/D^2)<1
%                 g(ii,jj)=g(ii,jj)+exp(2*c)*(exp(c*(D(kk)^2*(1/(x1(jj).^2-D(kk).^2))+E(kk)^2*(1/(y1(ii).^2-E(kk).^2)))));
%             end
%
%         end
%     end
% end