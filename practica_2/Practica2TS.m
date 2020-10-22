%% Teoria de Senyal: Pràctica 2 de laboratori:
%% Damià Casas & Pau Manyer

%% Estudi Previ:
clear all;
clc
close all;
% 2.6)
k0=3; N=20;
Omega=2*pi*(k0/N);

f=@(n) [cos(Omega*n)];
y=f([0:52]);

v=fft(y,N);

figure(1)
stem(y)
grid on

for k=0:N-1
    a=(1-exp(1i*2*pi*(k0-k)))/(1-exp(1i*2*(pi/N)*(k0-k)));
    b=(1-exp(-1i*2*pi*(k0+k)))/(1-exp(-1i*2*(pi/N)*(k0+k)));
    XN(k+1)=0.5*(a+b);
end

figure(2)
stem(v)
grid on
figure(3)
stem(XN)
grid on

error=abs(v-XN);
figure(4)
stem(error)
grid on

%%

clear all
clc

h=randn(1,1e3);
x=randn(1,1e7);
N=2.^[10:20];


%stem(y1);

y1=conv(x,h);
y2=process(x,h,N(3));

 error = abs(y1-y2);
 max(error)




%%

% Si X es una matriz, fft(X) trata las columnas de X como vectores y
% devuelve la transformada de Fourier de cada columna. Tornarà la FFT de
% tants punts com hi hagin a les columnes.

% If Y is a matrix, then ifft(Y) returns the inverse transform of each column of the matrix.

function yb=cc(xb,H) %Convolucio del bloc N=Lh+Lx-1
N=length(H);
Lx=length(xb);
xb=[xb zeros(1,N-Lx)];
X=fft(xb,N);
Y=X.*H;
yb=ifft(Y,N);
end

function y=process(x,h,N)
% INPUT: N: COM A MINIM HA DE SER length(x)+length(h)-1
L=length(h);
Lx=length(x);
M=N-L+1;  % Mmax
h=[h zeros(1,N-L)];
H=fft(h,N);  %llargada L+M-1

P=ceil(Lx/M);
longY = L+Lx-1;
y=zeros(1,longY);
for r=1:P-2 %Per cada bloc
    xi=[x([(r-1)*M+1:r*M])];
    yb=cc(xi,H); % cc ja s'encarregara de posar 0s al bloc i fer la convolucio
    A=(r-1)*M;
    
    temp = y([A+1:A+length(yb)]); %TODO: optimitzar per no calcular length
    temp = temp + yb;
    y([A+1:A+length(yb)]) = temp;
    disp(r)
end
xi=x([(P-1)*M+1:end]);
yb=cc(xi,H);  %llargada N
temp = y([(P-1)*M+1:end]);
temp = temp + yb(1:length(temp));
y([(P-1)*M+1:end]) = temp;
end






