function [tfd,A,am1] = myMDD(sig,c,th,disp,caseEEG)

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

N=length(sig);
am1 =wvd1(sig,N);
% th=0.3;

%% computation of theta
theta = -90:2:90;
[R,xp] = radon(abs(am1),theta);
pos=find(xp==0);
% R(pos,:)=R(pos,:)/sum(R(pos,:));
Rad_0=Dtrend(R(pos,:),10);
R_trend = R(pos,:)-Rad_0';
if disp==1
    figure;plot(90-theta,R(pos,:));xlabel('angle(degree)');ylabel('Magnitude');hold on;
    plot(90-theta,R_trend,'-r');xlabel('angle(degree)');ylabel('Magnitude');
end
Rad_0 = Rad_0/max(Rad_0);
% Rad_0=R(pos,:)/max(R(pos,:));
if disp==1
    figure;plot(90-theta,Rad_0);xlabel('angle(degree)');ylabel('Magnitude')
    hold on;plot(90-theta,th*ones(size(theta)),'-r');
end
[ peak,loc ] = Detection_pic(Rad_0);
pos_pic = find(peak>th);

if caseEEG
    theta = [90 0];
    E = 0.4*ones(size(theta));
else
    theta = 90-theta(loc(pos_pic));
    E = peak(pos_pic)/2;
end
% E=0.08*ones(size(theta));
% E=[0.07 0.08 0.15];
% E=0.27;
D = 0.01*ones(size(theta));
%% processing
A=kernel_MDD(sig,c,D,E, theta);
am=am1.*A;
tfd = (fft(ifftshift(am,1), [], 1));
tfd=  ifft(fftshift(tfd,2), [], 2);
tfd=real(tfd);
