clear all
close all

y = cos((2*pi/(4))*[0:20]);

figure()
stem(abs(fft(y)))
figure()
stem(abs(fft(y,40)))
