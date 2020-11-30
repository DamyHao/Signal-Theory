%% Teoria de Senyal: Pràctica 2 de laboratori:
%% Damià Casas & Pau Manyer

%% Estudi Previ:
clear all;
clc
close all;
% 2.6)
k0=3; N=32;
Omega=2*pi*(k0/N);
phi = pi/4;
f=@(n)(cos(Omega*n + phi));
y=f([0:64]);
% Recordar que k va de 0 a N-1
v=fft(y,N);

figure(1)
stem( y)
grid on
title("Senyal")
print('./eps/primerDibSenyal','-depsc','-tiff')
count = 1;
XN = zeros(1, N);

funcio = @(k)(0.5*((1-(exp(1i*2*(pi/N)*(k0-k)) + exp(-1i*2*(pi/N)*(k0+k))).^(N-1))./( 1- (exp(1i*2*(pi/N)*(k0-k)) + exp(-1i*2*(pi/N)*(k0+k))))));

for k=0:N
    %a=(1-exp(1i*2*pi*(k0-k)))/(1-exp(1i*2*(pi/N)*(k0-k)));
    %b=(1-exp(-1i*2*pi*(k0+k)))/(1-exp(-1i*2*(pi/N)*(k0+k)));
    XN(count)=cosDFT(k, k0, N, phi);
    count = count + 1;
end
ks=0:N;

figure(2)
stem([0:N-1] ,abs(v))
grid on
title("Amb fft de matlab")
ylabel('Valor absolut')
xlabel('k')
print('./eps/ambNostre','-depsc','-tiff')
figure(3)
stem([0:N], abs(XN))
grid on
title("Amb formula trobada")
ylabel('Valor absolut')
xlabel('k')
print('./eps/ambNostre','-depsc','-tiff')


%%

% 2.1)

clear all
close all
clc

h=randn(1,1e3);
x=randn(1,1e7);

expo=[10:20];
N=2.^expo;
for k=1:length(N)
    tic; y1=conv(x,h); tlist1(k)=toc;
    tic; y2=process(x,h,N(k)); tlist2(k)=toc;
    error = abs(y1-y2);
    maxerrorlist(k)=max(abs(error));
end

figure(5)
plot(expo,tlist1)
grid on
hold on
plot([10:20],tlist2)
xlabel("\nu")
ylabel("s")
title("Comparació de la durada de computació")
legend("conv()", "process()", "Location", "best")
print('./eps/durada12','-depsc','-tiff')


figure(6)
semilogy(expo,maxerrorlist)
grid on
xlabel("\nu")
ylabel("Max(abs(y_1-y_2))")
title("Màxim error absolut")
print('./eps/error','-depsc','-tiff')


%%
close all

L = 1000;
funcio23 = @(nu)((4*(2.^nu).*(1+nu))./(2.^nu - L + 1));
plot(expo, funcio23(expo))
grid on;
title("Multiplicacions necesaries per process");
legend("conv()", "Location", "best")
print -deps aprox23


%%
% 2.2)

clear all
clc
close all
% Fem el senyal 32
load('./DIAL_data/N32.MAT')

% Dibuixeu el senyal
figure(7)
stem(x)
grid on
xlabel("x")
ylabel("y")
title("Senyal")
print('./eps/dibuixSenyal','-dpdf', '-bestfit')


% Quina és la duració del senyal:
fs=8e3; Ts=1/fs;
L=length(x);
N = L;
display(L,'la llargada de x és')
display(L,'La duració del senyal en número de mostres és')
display(Ts*L,'La duració del senyal en segons és')

% Escoltem el senyal:
player = audioplayer(x,fs);
play(player);


% DFT i dibuix del mòdul:

X=fft(x);
figure(8)
plot(abs(X))
grid on
print('./eps/dft','-depsc', '-tiff')

% Trobi els índexs (valors de k) per a cadascun dels 4 pics de la DFT.
Mod=abs(X);
[values,loc]=sort(Mod);
klist=[values(end-3) values(end-2) values(end-1) values(end)];
locationPics=[loc(end-3) loc(end-2) loc(end-1) loc(end)];

% Obtingui les freqüències digitals corresponents als dos pics detectats
Omega1=2*(pi/N)*locationPics(1);
Omega2=2*(pi/N)*locationPics(3);

discreteOmegas= 2*(pi/N).*locationPics;
% Quines són les freqüències analògiques (en Hz) corresponents a les freqüències
% digitals anteriors?
analogOmegas = discreteOmegas.*fs;
analogFreqs = analogOmegas./(2*pi);

w1=Omega1*fs;
w2=Omega2*fs;

f1=w1/(2*pi);
f2=w2/(2*pi);

% Quina és la tecla corresponent al senyal x utilitzant la Figura 1?

%Tallem la senyal:
onTallar1 = 3/4*L;
onTallar2 = 1/2*L;

senyal1 = x([1:onTallar1]);
senyal2 = x([1:onTallar2]);

dft1 = fft(senyal1, N);
dft2 = fft(senyal2, N);

figure()
plot(abs(dft1))
xlabel('k')
title('DFT de 3/4 del senyal')
grid on
print('./eps/dft1','-depsc', '-tiff')

figure()
plot(abs(dft2))
xlabel('k')
title('DFT de 1/2 del senyal')
grid on
print('./eps/dft2','-depsc', '-tiff')

zeroPad=fft(x, 2^16);
figure()
stem(abs(zeroPad), 'Marker', 'none')
grid on

%% Ultim apartat 2.3.3
close all;
N = 1024;
n = [1:N];
hamm = hamming(N);
hammingWindowing = hamm'.*x;

figure()
stem(hamm);
xlabel('n');
ylabel('w[n]')
title('Hamming windowing mida N');
print('./eps/hamming','-depsc', '-tiff')

figure();
stem(hammingWindowing);
xlabel('n');
ylabel('x[n]')
title('Hamming windowing mida N aplicat a la senyal.');
print('./eps/hammingSignal','-depsc', '-tiff')

frecuencyDomainBartlett = fft(hamm);
figure()
stem(20*abs(frecuencyDomainBartlett), 'Marker', 'none')
set(gca,'yscal','log')
%semilogy(n, abs(frecuencyDomainBartlett))
ylabel('dB');
xlabel('k');
title('Fft Hamming');
print('./eps/hammingFreq','-depsc', '-tiff')
% Per veure-ho millor:
frecuencyDomainBartlett = fft(hamm, 2^16);
figure()
%stem(20*abs(frecuencyDomainBartlett), 'Marker', 'none')
%set(gca,'yscal','log')
semilogy([1:1:2^16], 20*abs(frecuencyDomainBartlett))
ylabel('dB');
xlabel('k');
title('Fft Hamming');
print('./eps/hammingFreqP','-depsc', '-tiff')

barletFft = fft(hammingWindowing);
figure;
stem(abs(barletFft), 'Marker', 'none');
xlabel('k');
title('FFT de windowed signal amb Hamming de N');
print('./eps/hammingSignalFFt','-depsc', '-tiff');

% Ara amb bartlett
hamm = bartlett(N);
hammingWindowing = hamm'.*x;

figure()
stem(hamm);
xlabel('n');
ylabel('w[n]')
title('Bartlett windowing mida N');
print('./eps/bartlett','-depsc', '-tiff')

figure();
stem(hammingWindowing);
xlabel('n');
ylabel('x[n]')
title('Bartlett windowing mida N aplicat a la senyal.');
print('./eps/bartlettSignal','-depsc', '-tiff')

frecuencyDomainBartlett = fft(hamm);
figure()
stem(20*abs(frecuencyDomainBartlett), 'Marker', 'none')
set(gca,'yscal','log')
%semilogy(n, abs(frecuencyDomainBartlett))
ylabel('dB');
xlabel('k');
title('Fft Bartlett');
print('./eps/bartlettFreq','-depsc', '-tiff')
% Per veure-ho millor:
frecuencyDomainBartlett = fft(hamm, 2^16);
figure()
%stem(20*abs(frecuencyDomainBartlett), 'Marker', 'none')
%set(gca,'yscal','log')
semilogy([1:1:2^16], 20*abs(frecuencyDomainBartlett))
ylabel('dB');
xlabel('k');
title('Fft bartlett');
print('./eps/bartlettFreqP','-depsc', '-tiff')

barletFft = fft(hammingWindowing);
figure;
stem(abs(barletFft), 'Marker', 'none');
xlabel('k');
title('FFT de windowed signal amb Bartlett de N');
print('./eps/bartlettSignalFFt','-depsc', '-tiff');


%Ara partint la senyal fins a 1/4:

N = 256;
XX = x(1:256);
n = [1:N];
hamm = hamming(N);
hammingWindowing = hamm'.*XX;

figure()
stem(hamm);
xlabel('n');
ylabel('w[n]')
title('Hamming windowing mida N');

figure();
stem(hammingWindowing);
xlabel('n');
ylabel('x[n]')
title('Hamming windowing mida N aplicat a la senyal.');

frecuencyDomainBartlett = fft(hamm);
figure()
stem(20*abs(frecuencyDomainBartlett), 'Marker', 'none')
set(gca,'yscal','log')
%semilogy(n, abs(frecuencyDomainBartlett))
ylabel('dB');
xlabel('k');
title('Fft Hamming');
% Per veure-ho millor:
frecuencyDomainBartlett = fft(hamm, 2^16);
figure()
%stem(20*abs(frecuencyDomainBartlett), 'Marker', 'none')
%set(gca,'yscal','log')
semilogy([1:1:2^16], 20*abs(frecuencyDomainBartlett))
ylabel('dB');
xlabel('k');
title('Fft Hamming');

barletFft = fft(hammingWindowing);
figure;
stem(abs(barletFft), 'Marker', 'none');
xlabel('k');
title('FFT de windowed signal amb Hamming de N');
print('./eps/hammingSignalFFtW','-depsc', '-tiff');

% Ara amb bartlett
hamm = bartlett(N);
hammingWindowing = hamm'.*XX;

figure()
stem(hamm);
xlabel('n');
ylabel('w[n]')
title('Bartlett windowing mida N');

figure();
stem(hammingWindowing);
xlabel('n');
ylabel('x[n]')
title('Bartlett windowing mida N aplicat a la senyal.');

frecuencyDomainBartlett = fft(hamm);
figure()
stem(20*abs(frecuencyDomainBartlett), 'Marker', 'none')
set(gca,'yscal','log')
%semilogy(n, abs(frecuencyDomainBartlett))
ylabel('dB');
xlabel('k');
title('Fft Bartlett');
% Per veure-ho millor:
frecuencyDomainBartlett = fft(hamm, 2^16);
figure()
%stem(20*abs(frecuencyDomainBartlett), 'Marker', 'none')
%set(gca,'yscal','log')
semilogy([1:1:2^16], 20*abs(frecuencyDomainBartlett))
ylabel('dB');
xlabel('k');
title('Fft bartlett');

barletFft = fft(hammingWindowing);
figure;
stem(abs(barletFft), 'Marker', 'none');
xlabel('k');
title('FFT de windowed signal amb Bartlett de N');
print('./eps/bartlettSignalFFtWW','-depsc', '-tiff');

%%

% Si X es una matriz, fft(X) trata las columnas de X como vectores y
% devuelve la transformada de Fourier de cada columna. Tornarà la FFT de
% tants punts com hi hagin a les columnes.

% If Y is a matrix, then ifft(Y) returns the inverse transform of each column of the matrix.

function yb=cc(xb,H) 
% Convolucio circular d'ordre N=Lh+Lx-1 (implicita en llargada de H). Helper de process.
% INPUT:
%   xb: bloc de mida M.
%   H: Fourier transform de respota impulsional h(n). Ha de ser allargat
%       amb 0s fins N
% OUTPUT:
%   yb: convolucio circular

N=length(H);
xb=[xb zeros(1,N-length(xb))]; % Fem el vector més llarg
X=fft(xb,N); % Fast Fourier trransform
Y=X.*H;
yb=ifft(Y,N); %La inversa de la fft tornara el producte circular. Si 
end

function y=process(x,h,N)
% Calcula la senyal de sortida de un sistema FIR trencant-la en blocs i
% fent la fast fourirer transform.
% INPUT:
%   x: senyal de entrada.
%   h: resposta impulsional
%   N: Ordre de la FFT.
% OUTPUT:
%   y: senyal de sortida
L=length(h);   % Calculem la longitud del vector h
Lx=length(x);
M=N-L+1;  % Calculem la maxima longitud de M. Si fos mes petit, 
% al fer la inversa dona la convolucio circular i no seria igual a la lineal
h=[h zeros(1,N-L)];
%H=fft(h,N);
H=fft(h,N);  % Calculem la DFT de N punts de h, tindra llargada L+M-1
P=ceil(Lx/M);  % Calculem P i ens assegurem que es enter arrodonint cap a dalt
% El ultim bloc sera mes petit que M
longY = L+Lx-1; % llargada que tindra la convolucio final
y=zeros(1,longY);
for r=1:P-1 % Per els P-1 primers blocs de llargada M
    xi=x((r-1)*M+1:r*M);
    yi=cc(xi,H); % cc calcula la convolucio (circular = lineal) entre el bloc xi i H
    A=(r-1)*M;
    temp = y(A+1:A+length(yi));    % Aquest troç concatena adequadament els
    temp = temp + yi;              % blocs yi
    y([A+1:A+length(yi)]) = temp;
end
% Ultim bloc:
xi=x([(P-1)*M+1:end]);  % per els valors sobrants de x tambe calculem la
yi=cc(xi,H);            % convolucio amb H i concatenem y
temp = y((P-1)*M+1:end);
temp = temp + yi(1:length(temp));
y((P-1)*M+1:end) = temp;

end


function y=cosDFT(k,k0, N, phi)
    if (k == k0) 
        y = (N/2)*exp(1i*phi);
    elseif (k == (N - k0))
        y = (N/2)*exp(-1i*phi);
    else
        y = 0;
    end
end


function v = berlettWindowing(n,N)
    v = (N-1)/2 - abs(n-(N-1)./2);
end

