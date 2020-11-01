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

% 2.1)

clear all
clc

h=randn(1,1e3);
x=randn(1,1e7);

%stem(y1);

y1=conv(x,h);
y2=process(x,h,2^(15));

expo=[11:20];
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
plot([11:20],tlist2)
xlabel("\nu")
ylabel("Durada de computació")
title("Comparació de la durada de computació de cada mètode en funció de N (=2^\nu)")

figure(6)
semilogy(expo,maxerrorlist)
grid on
xlabel("\nu")
ylabel("Max(abs(y_1-y_2))")
title("Màxim error abs(y_1-y_2) en funció de l'exponent de N (=2^\nu)")

%% Finestra rectangular de llargada N:

clear all
clc

format long g

% Transformada de S:
S=@(omega) [pi*(dirac(omega-2*pi/3)+dirac(omega-4*pi/3))];
omega=linspace(0,2*pi,961);

% Representem les deltas de Dirac amb valor 1 per comoditat.
TransfS=S(omega);
[loc,index]=sort(TransfS);
TransfS(index(end)) = 1; TransfS(index(end-1)) = 1;     

figure(1)
plot(omega,TransfS,'-o')
grid on
xlabel('Omega')
ylabel('Valor transformada')
title('Transformada de Fourier de s[n]')
axis([-1 7.5 -0.2 1.2])

% Transformada de W:
N=20;
W=@(omega) [((sin(omega.*N/2))./(sin(omega/2))).*exp(-1i*omega*(N-1)/2)];
omega=linspace(0,2*pi,961);   

TransfW=W(omega);
TransfW=TransfW(1:end-1);

figure(2)
plot(omega(1:end-1),abs(TransfW))
grid on
xlabel('Omega')
ylabel('Valor transformada')
title('Transformada de Fourier de w[n]')
axis([0 2*pi -0.2 max(abs(TransfW))+3])

% Transformada de X:
X=@(omega) [((sin((omega-2*pi/3).*N/2))./(sin((omega-...
    2*pi/3)/2))).*exp(-1i*(omega-2*pi/3)*(N-1)/2)+...
    ((sin((omega-4*pi/3).*N/2))./(sin((omega-4*pi/3)/2))).*exp(...
    -1i*(omega-4*pi/3)*(N-1)/2)];

omega=linspace(0,2*pi,961);

TransfX=X(omega);
TransfX=TransfX(1:end-1);

figure(3)
plot(omega(1:end-1),abs(TransfX))
grid on
xlabel('Omega')
ylabel('Valor transformada')
title('Transformada de Fourier de x[n]')
axis([0 2*pi -0.2 max(abs(TransfX))+3])

% DFT:
for k=0:N-1
    DFT(k+1)=abs(X(0.000001+(2*pi*k/N)));
end

figure(4)
stem(DFT)
grid on
xlabel('k')
ylabel('Valor DFT')
title('DFT de x[n]')

%% Finestra rectangular fixa (zero padding):

clear all
clc

format long g

% Transformada de S:
S=@(omega) [pi*(dirac(omega-2*pi/3)+dirac(omega-4*pi/3))];
omega=linspace(0,2*pi,961);

% Representem les deltas de Dirac amb valor 1 per comoditat.
TransfS=S(omega);
[loc,index]=sort(TransfS);
TransfS(index(end)) = 1; TransfS(index(end-1)) = 1;     

figure(1)
plot(omega,TransfS,'-o')
grid on
xlabel('Omega')
ylabel('Valor transformada')
title('Transformada de Fourier de s[n]')
axis([-1 7.5 -0.2 1.2])

% Transformada de W:
Nfinestra=20;
W=@(omega) [((sin(omega.*Nfinestra/2))./(sin(omega/2))).*exp(-1i*omega*(Nfinestra-1)/2)];
omega=linspace(0,2*pi,961);   

TransfW=W(omega);
TransfW=TransfW(1:end-1);

figure(2)
plot(omega(1:end-1),abs(TransfW))
grid on
xlabel('Omega')
ylabel('Valor transformada')
title('Transformada de Fourier de w[n]')
axis([0 2*pi -0.2 max(abs(TransfW))+3])

% Transformada de X:
X=@(omega) [((sin((omega-2*pi/3).*Nfinestra/2))./(sin((omega-...
    2*pi/3)/2))).*exp(-1i*(omega-2*pi/3)*(Nfinestra-1)/2)+...
    ((sin((omega-4*pi/3).*Nfinestra/2))./(sin((omega-4*pi/3)/2))).*exp(...
    -1i*(omega-4*pi/3)*(Nfinestra-1)/2)];

omega=linspace(0,2*pi,961);

TransfX=X(omega);
TransfX=TransfX(1:end-1);

figure(3)
plot(omega(1:end-1),abs(TransfX))
grid on
xlabel('Omega')
ylabel('Valor transformada')
title('Transformada de Fourier de x[n]')
axis([0 2*pi -0.2 max(abs(TransfX))+3])

% DFT:
N=20;
for k=0:N-1
    DFT(k+1)=abs(X(0.000001+(2*pi*k/N)));
end

figure(4)
stem(DFT)
grid on
xlabel('k')
ylabel('Valor DFT')
title('DFT de x[n]')

%% Finestra de Hann

clear all
clc

format long g

% Transformada de S:
S=@(omega) [pi*(dirac(omega-2*pi/3)+dirac(omega-4*pi/3))];
omega=linspace(0,2*pi,961);

% Representem les deltas de Dirac amb valor 1 per comoditat.
TransfS=S(omega);
[loc,index]=sort(TransfS);
TransfS(index(end)) = 1; TransfS(index(end-1)) = 1;     

figure(1)
plot(omega,TransfS,'-o')
grid on
xlabel('Omega')
ylabel('Valor transformada')
title('Transformada de Fourier de s[n]')
axis([-1 7.5 -0.2 1.2])

% Transformada de W:
N=20;
W=@(oo)(pi.*((sin(oo.*N/2))./(sin(oo./2))).*...
    exp(-1i*oo.*(N-1)./2)-pi.*((sin((oo-2*pi/N).*N/2))./...
    (sin((oo-2*pi/N)./2))).*exp(-1i*(oo-2*pi/N).*(N-1)/2)...
    -pi.*((sin((oo+2*pi./N).*N./2))./(sin((oo+2*pi./N)./2))).*...
    exp(-1i*(oo+2*pi./N).*(N-1)./2));

omega=linspace(0,2*pi,961);   

TransfW=W(omega);
TransfW=TransfW(1:end-1);

figure(2)
plot(omega(1:end-1),abs(TransfW))
grid on
xlabel('Omega')
ylabel('Valor transformada')
title('Transformada de Fourier de w[n]')
axis([0 2*pi -0.2 max(abs(TransfW))+3])

% Transformada de X:
X=@(omega) [((sin((omega-2*pi/3).*N/2))./(sin((omega-...
    2*pi/3)/2))).*exp(-1i*(omega-2*pi/3)*(N-1)/2)+...
    ((sin((omega-4*pi/3).*N/2))./(sin((omega-4*pi/3)/2))).*exp(...
    -1i*(omega-4*pi/3)*(N-1)/2)];

omega=linspace(0,2*pi,961);

TransfX=X(omega);
TransfX=TransfX(1:end-1);

figure(3)
plot(omega(1:end-1),abs(TransfX))
grid on
xlabel('Omega')
ylabel('Valor transformada')
title('Transformada de Fourier de x[n]')
axis([0 2*pi -0.2 max(abs(TransfX))+3])

% DFT:
N=20;
for k=0:N-1
    DFT(k+1)=abs(X(0.000001+(2*pi*k/N)));
end

figure(4)
stem(DFT)
grid on
xlabel('k')
ylabel('Valor DFT')
title('DFT de x[n]')







%%
% 2.2)

clear all
clc
load('N32.MAT')

% Dibuixeu el senyal
figure(7)
stem(x)
grid on 

figure(8)
stem(x(1:100))
grid on 

% Quina és la duració del senyal:
fs=8e3; Ts=1/fs;
L=length(x);
display(L,'la llargada de x és')
display(L,'La duració del senyal en número de mostres és')
display(Ts*L,'La duració del senyal en segons de mostres és')

% DFT i dibuix del mòdul:
N=2^8;
X=fft(x,N);
figure(9)
plot(abs(X))
grid on
xlabel('k')
ylabel('X[k]')
title('DFT[x]=X[k]')

% Trobi els índexs (valors de k) per a cadascun dels 4 pics de la DFT.
Mod=abs(X);
[loc,index]=sort(Mod);
klist=[index(end-3) index(end-2) index(end-1) index(end)];
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
L=length(h);   % Calculem la longitud del vector h
Lx=length(x);
M=N-L+1;  % Calculem la maxima longitud de M
h=[h zeros(1,N-L)];
H=fft(h,N);  % Calculem la DFT de N punts de h, tindra llargada L+M-1
P=ceil(Lx/M);  % Calculem P i ens assegurem que es enter amb el comando ceil()
longY = L+Lx-1;
y=zeros(1,longY);
for r=1:P-1         % Per els P-1 primers blocs de llargada M fem:
    xi=x((r-1)*M+1:r*M);
    yi=cc(xi,H); % cc calcula la convolucio entre el bloc xi i H
    A=(r-1)*M; 
    temp = y(A+1:A+length(yi));    % Aquest troç concatena adequadament els 
    temp = temp + yi;              % blocs yi
    y([A+1:A+length(yi)]) = temp;
end
xi=x([(P-1)*M+1:end]);  % per els valors sobrants de x tambe calculem la  
yi=cc(xi,H);            % convolucio amb H i concatenem y
temp = y((P-1)*M+1:end); 
temp = temp + yi(1:length(temp));
y((P-1)*M+1:end) = temp;
end






