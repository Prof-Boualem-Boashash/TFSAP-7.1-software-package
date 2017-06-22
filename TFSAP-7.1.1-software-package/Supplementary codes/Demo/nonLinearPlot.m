function nonLinearPlot(tfr,f,t)

N = size(tfr,2);

figure;imagesc(1:N,t,tfr);axis xy;colormap((jet))

ylabel('Time(sec)');xlabel('Frequency(Hz)'); title('');


YtickLabel = get(gca,'XtickLabel');
for k=1:length(YtickLabel)
    YtickLabel2(k) = str2num(YtickLabel{k});
end

YtickLabel_new = 0.01*round(100*(f(YtickLabel2)));

YtickLabel_new = num2str(YtickLabel_new);

set(gca,'XtickLabel',YtickLabel_new);

end