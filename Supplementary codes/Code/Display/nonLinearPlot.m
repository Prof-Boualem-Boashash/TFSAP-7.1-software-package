function nonLinearPlot(tfr,f,t)

N = size(tfr,1);

figure;imagesc(t,1:N,tfr);axis xy;colorbar;colormap((jet))

xlabel('Time(sec)');ylabel('Frequency(Hz)'); title('');


YtickLabel = get(gca,'YtickLabel');
for k=1:length(YtickLabel)
    YtickLabel2(k) = str2num(YtickLabel{k});
end

YtickLabel_new = (f(YtickLabel2));

YtickLabel_new = num2str(YtickLabel_new);

set(gca,'YtickLabel',YtickLabel_new);

end