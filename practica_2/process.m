function y=process(x,h,N)
% INPUT: N: COM A MINIM HA DE SER length(x)+length(h)-1
% Calcula la senyal de sortida de un sistema FIR trencant-la en blocs i
% fent la fast fourirer transform
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
