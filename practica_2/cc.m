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