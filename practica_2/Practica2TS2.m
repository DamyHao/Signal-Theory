%% Teoria de Senyal: Pràctica 2 de laboratori:
%% Damià Casas & Pau Manyer

%% Estudi Previ:
clear all;
clc
close all;
% 2.6)
k0=3; N=16;
Omega=2*pi*(k0/N);

f=@(n) [cos(Omega*n)];
y=f([0:52]);

v=fft(y,N);

figure(1)
stem(y)
grid on
title("Senyal")
count = 1;
XN = zeros(1, N);

funcio = @(k)(0.5*((1-(exp(1i*2*(pi/N)*(k0-k)) + exp(-1i*2*(pi/N)*(k0+k))).^(N-1))./( 1- (exp(1i*2*(pi/N)*(k0-k)) + exp(-1i*2*(pi/N)*(k0+k))))));

for k=0:N
    a=(1-exp(1i*2*pi*(k0-k)))/(1-exp(1i*2*(pi/N)*(k0-k)));
    b=(1-exp(-1i*2*pi*(k0+k)))/(1-exp(-1i*2*(pi/N)*(k0+k)));
    XN(count)=0.5*(a+b);
    count = count + 1;
end
ks=0:N;
XN = funcio(ks);

figure(2)
stem(abs(v))
grid on
title("Amb fft")
figure(3)
stem(abs(XN))
grid on
title("Amb formula trobada")


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
ylabel("Durada de computació en segons")
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
% Fem el senyal 32
load('./DIAL_data/N32.MAT')

% Dibuixeu el senyal
figure(7)
stem(x)
grid on
xlabel("x")
ylabel("y")
title("Senyal")
print('./eps/dibuixSenyal','-depsc','-tiff')


% Quina és la duració del senyal:
fs=8e3; Ts=1/fs;
L=length(x);
display(L,'la llargada de x és')
display(L,'La duració del senyal en número de mostres és')
display(Ts*L,'La duració del senyal en segons és')

% Escoltem el senyal:
player = audioplayer(x,fs);
play(player);


% DFT i dibuix del mòdul:
N=2^8;
X=fft(x,N);
figure(8)
plot(abs(X))
grid on
print('./eps/dft','-depsc','-tiff')


% Trobi els índexs (valors de k) per a cadascun dels 4 pics de la DFT.
Mod=abs(X);
[values,loc]=sort(Mod);
klist=[values(end-3) values(end-2) values(end-1) values(end)];
valorspics=[loc(end-3) loc(end-2) loc(end-1) loc(end)];

% Obtingui les freqüències digitals corresponents als dos pics detectats
Omega1=2*(pi/N)*klist(1);
Omega2=2*(pi/N)*klist(3);

% Quines són les freqüències analògiques (en Hz) corresponents a les freqüències
% digitals anteriors?

w1=Omega1*fs;
w2=Omega2*fs;

f1=w1/(2*pi);
f2=w2/(2*pi);

% Quina és la tecla corresponent al senyal x utilitzant la Figura 1?



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





