load La4;
soundsc(Nota,Fs);
figure(1)
NFFT=2^ceil(log2(length(Nota)));
TF=fft(Nota,NFFT);
 
f=[0:1/NFFT:0.5]*Fs;
plot(f,abs(TF(1:NFFT/2+1)),'r')
print('./eps/original','-depsc', '-tiff');
axis([0 8000 get(gca,'YLim')])
M=4;
xD_Filter=decimate(Nota,M,'FIR');

%%
soundsc(xD_Filter,Fs/M);
%pause
figure(2)
TF1=fft(xD_Filter,NFFT);
plot(f/M,abs(TF1(1:NFFT/2+1)),'r')
axis([0 8000 get(gca,'YLim')])

print('./eps/ambFiltreHamming','-depsc', '-tiff');



%%
M=4;
xD_NoFilter=downsample(Nota,M);
soundsc(xD_NoFilter,Fs/M);
%pause
figure(3)
TF=fft(xD_NoFilter,NFFT);
plot(f/M,abs(TF(1:NFFT/2+1)),'r')
axis([0 8000 get(gca,'YLim')])
print('./eps/senseFiltre','-depsc', '-tiff');


figure(4)
plot(f/M,abs(TF(1:NFFT/2+1)),'r')
hold on;
plot(f/M,abs(TF1(1:NFFT/2+1)),'b')
hold off;
axis([0 8000 get(gca,'YLim')])
print('./eps/sobreposicio','-depsc', '-tiff');

