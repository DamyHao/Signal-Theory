clear all
close all

% y = cos((2*pi/(4))*[0:20]);
% 
% figure()
% stem(abs(fft(y)))
% figure()
% stem(abs(fft(y,40)))
s=tf('s');
options = bodeoptions;
options.FreqUnits = 'Hz';
RC = 1e-4;
H=1/(RC*s + 1);
syms s
RC = 100*1e-9;
T=1/(RC*s + 1);
bode(H, options);grid on %Bode teoric
hold on;
%%
close all
x = [1 9 13 14 15 15.42 16 17 100 10000];
x = x.*100*2*pi;
volt = [1010 864 775 754 731 725 706 633 150 3.96];
volt= volt./1000;

figure()
semilogx(x, volt, 'o');
hold on;
RC = 100*1e-9*1000;
modulH = @(w)(1./(sqrt((RC.*w).^2+1)));

%CANVIAR PER LOGLOG
semilogx(10.^[1:0.2:8], modulH(10.^[1:0.2:8]))
hold off;
legend("Resposta experimental", "Resposta teorica");
ylabel("Modulus")
xlabel("\omega");

figure()
angles = [8 24 41 43 44 45 46 58 82 90];
angles = angles*(pi/180)
%CANVIAR PER LOGLOG
semilogx(x, -angles, 'o');
hold on;
hAngle = @(s)(angle(1./(1i.*RC.*s + 1)));
semilogx(10.^[1:0.2:8], hAngle(10.^[1:0.2:8]))
hold off;
legend("Resposta experimental", "Resposta teorica");
xlabel("\omega");
ylabel("argument");
