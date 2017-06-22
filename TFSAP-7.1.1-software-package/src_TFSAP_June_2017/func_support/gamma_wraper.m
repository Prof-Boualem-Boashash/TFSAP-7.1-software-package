function [f_filter] = gamma_wraper ( pars )

% This function is a wrapper for calling gamma-complex function  

beta = pars(1);
window_length = pars(2);

v=-.5:1/(window_length-1):.5;
s=beta+(1i*pi*v);
f = (abs(gammacomplex(s)).^2)/abs(gammacomplex(beta)).^2;

f = abs(f);

if mod(length(f),2)
    N = floor(length(f)/2);
else
    N = length(f)/2;
end

f_filter=f(N+1:end);
f_filter=[f_filter f(1:N)];