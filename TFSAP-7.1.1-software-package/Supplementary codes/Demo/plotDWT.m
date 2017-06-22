function plotDWT(sig,wavname)

[LO_R,HI_R] = mywfilters(wavname);
[signal_dec, dimplan] = mymallat1D(LO_R,HI_R,sig);
fs = 32;
n=1;
deb = 1;
t=0:1/(fs):(length(sig)-1)/(fs);
% figure;plot(t,sig);xlabel('Time (s)');ylabel('Magnitude')
signal_dim=signal_dec(deb:deb+dimplan(n)-1);
deb = deb+dimplan(n);
t1=0:1/(fs/2):(length(signal_dim)-1)/(fs/2);
% figure;plot(t1,signal_dim);xlabel('Time (s)');ylabel('Magnitude')
n=2;
signal_dim2=signal_dec(deb:deb+dimplan(n)-1);
deb = deb+dimplan(n);
t2=0:1/(fs/4):(length(signal_dim2)-1)/(fs/4);
% figure;plot(t2,signal_dim2);xlabel('Time (s)');ylabel('Magnitude')
n=3;
signal_dim3=signal_dec(deb:deb+dimplan(n)-1);
deb = deb+dimplan(n);
t3=0:1/(fs/8):(length(signal_dim3)-1)/(fs/8);
% figure;plot(t3,signal_dim3);xlabel('Time (s)');ylabel('Magnitude')

figure
subplot(2,2,1)
plot(t,sig);axis([0 8 min(sig) max(sig)]);xlabel('Time (s)')
title('Original signal')

subplot(2,2,2)
plot(t1,signal_dim);axis([0 8 min(signal_dim) max(signal_dim)]);xlabel('Time (s)')
title('Level 1 of the DWT decomposition')

subplot(2,2,3)
plot(t2,signal_dim2);axis([0 8 min(signal_dim2) max(signal_dim2)]);xlabel('Time (s)')
title('Level 2 of the DWT decomposition')

subplot(2,2,4)
plot(t3,signal_dim3);axis([0 8 min(signal_dim3) max(signal_dim3)]);xlabel('Time (s)')
title('Level 3 of the DWT decomposition')