%% Teoria de Senyal: Pràctica 1 de laboratori:
%% Damià Casas & Pau Manyer

%% Estudi Previ:
%
% Qüestió 1.3:
%
% $$y[n+2]=x[n+2]-\frac{2}{3}x[n+1]+\frac{7}{8}y[n+1]-\frac{3}{32}y[n] $$
%
% L'equació en diferències finites ens demana que el vector $y$ contingui
% almenys $y[0]$ i $y[1]$ per poder inicialitzar la recurrència.

close all
N = 10;
for k = 0:N-1
    x(k+1)=delta(k)-2*delta(k-1)+4*delta(k-2);
end
a=[1 ; -7/8 ; 3/32];
b=[1; -2/3];
N=length(x);

%
% function y=iir_fd1(a,b,x,N)
% y(1)=b(1)*x(1);
% y(2)=b(1)*x(2)+b(2)*x(1)-a(2)*y(1);
% for k=1:N-2
%     compox=b(1)*x(k+2)+b(2)*x(k+1);
%     compoy=-a(2)*y(k+1)-a(3)*y(k);
%     y(k+2)=(1/(a(1)))*(compox+compoy);
% end
% end

% function y=iir_fd1(a,b,x,N)
% y=[0 0];
% for k=3:N
%     compox=b(1)*x(k)+b(2)*x(k-1);
%     compoy=-a(2)*y(k-1)-a(3)*y(k-2);
%     y(k)=(1/(a(1)))*(compox+compoy);
% end
% end
%
% function y=iir_fd2(a,b,x,N)
% y(1)=b(1)*x(1);
% y(2)=b(1)*x(2)+b(2)*x(1)-a(2)*y(1);
% for k=1:N-2
%     y(k+2)=(1/(a(1)))*(b(1)*x(k+2)+b(2)*x(k+1)-a(2)*y(k+1)-a(3)*y(k));
% end
% end

% function y=iir_fd2(a,b,x,N)
% y=[0 0];
% for k=3:N
%     y(k)=(1/(a(1)))*(b(1)*x(k)+b(2)*x(k-1)-a(2)*y(k-1)-a(3)*y(k-2));
% end
% end

%% Activitats al laboratori:

%%
% Activitat 3.1:
%
% Generem el senyal $x[n]$ tal que
%
% $$x[n]=delta[n]-2delta[n-1]+4delta[n-2]$$
%

N=10; %Es va dir que es canvies N=40 per N=10
for k=0:N-1
    x(k+1)=delta(k)-2*delta(k-1)+4*delta(k-2);
end

%%
% Ara representem el senyal $x[n]$ després de passar-lo per les diferents
% funcions implementades $iir\_df1$ i $iir\_fd2$.
%
% El senyal resultant per $iir\_df1$ és el següent:

a=[1 -7/8 3/32];
b=[1 -2/3];
N=length(x);

y1=iir_fd1(a,b,x,N);

figure(1)
stem(y1,'-o')
grid on
title('iir\_df1(x[n])');

%%
% De la mateixa manera, representem el senyal resultant per $iir\_df1$:

y2=iir_fd2(a,b,x,N);

figure(2)
stem(y2,'-o')
grid on
title('iir\_df2(x[n])');

% Veiem per tant que el resultat es el mateix amb els dos mètodes

%%
% Activitat 3.2:
%
% El objectiu ara és el d'obtenir la resposta impulsional del filtre que
% hem implementat en les funcions $iir\_df1$ i $iir\_fd2$. Es a dir el
% filtre descrit per:
%
% $$y[n+2]=x[n+2]-\frac{2}{3}x[n+1]+\frac{7}{8}y[n+1]-\frac{3}{32}y[n] $$
%
% Per obtenir la respota impulsional $h_{IIR}[n]$, només hem de fer passar
% per el filtre un senyal equivalent a una delta tal que $x[n]=delta[n]$.

clear all
clc

%Delta de N
x=[1 zeros(1,39)];

a=[1 -7/8 3/32];
b=[1 -2/3];
N=length(x);

h=iir_fd1(a,b,x,N);

figure(3)
stem(h,'-o')
grid on
title('Resposta impulsional h_{IRR}[n] del filtre (1)');

%%
% Mirem ara si per un senyal graó obtenim la mateixa resposta quan
% el fem passar pel filtre IIR(1) (calcular $iir\_df1(x[n])$) i quan el fem
% passar per un sistema FIR amb la mateixa resposta impulsional que el
% anterior $h_{IIR}[n]$. Per el sistema FIR calcularem $x*h_{IIR}[n]$ on 
% * representa el producte de convolució Podem dir que $h_{IIR}[n]$ pertany 
% a un sistema FIR perquè es finita

clear x % treu la variable de la memoria
N=10;
% for k=0:N-1
%     x(k+1)=step(k);
% end
%x1 = zeros(1,40);
%x2 = zeros(1, 49);

for k=1:N
    x(k)=step(k);
end
for k=1:21
    x2(k)=step(k);
end

a=[1 ; -7/8 ; 3/32]; b=[1; -2/3];
h=h([1:10]);

y1=iir_fd1(a, b, x2, 21);
y2=conv(x, h); %Convolucio. El resultat sera un vector de llargada x.len + h.len - 1

figure(4)
stem(y1,'-b')
hold on
stem(y2,'-*')
grid on
title('');
legend('Resposta del filtre (1) IIR','Resposta del sistema FIR','Location','southwest');
hold off

%%
% 


%%
% Activitat 3.3:
%
clear all

N=40;

b=[1 -sqrt(2) 1];
a=[1 1.1 0.5];

x1=cos((pi/4)*[0:N]);
x2=cos((pi*3/4)*[0:N]);

% Filter: a(1) ha de ser 1, si no normalitzarà
% a() son els que multipliquen la sortida y
% b() son els que multipliquen la entrada x (veure apunts de filtres
% discrets)
y1=filter(b,a,x1); 
y2=filter(b,a,x2);

figure(5)
stem(y1,'-r')
hold on
stem(y2,'-b')
grid on
title('Sortida del filtre amb entrades cosenoidals');
legend('pi/4','3*pi/4');
hold off


%%
% Observem que atenua l'entrada del cosinus amb freqüència menor i deixa
% igual la de freqüencia major. Per tant en un principi apostariem per un
% filtre pas alt pero tambe es podria tractar de un pas banda.

%%

function y=iir_fd1(a,b,x,N)
y(1)=b(1)*x(1);
y(2)=b(1)*x(2)+b(2)*x(1)-a(2)*y(1);
for k=1:N-2
    compox=b(1)*x(k+2)+b(2)*x(k+1);
    compoy=-a(2)*y(k+1)-a(3)*y(k);
    y(k+2)=(1/(a(1)))*(compox+compoy);
end
end

function y=iir_fd2(a,b,x,N)
y(1)=b(1)*x(1);
y(2)=b(1)*x(2)+b(2)*x(1)-a(2)*y(1);
for k=1:N-2
    y(k+2)=(1/(a(1)))*(b(1)*x(k+2)+b(2)*x(k+1)-a(2)*y(k+1)-a(3)*y(k));
end
end

function out=delta(in)
if in==0
    out=1;
else
    out=0;
end
end

function out=step(in)
if in>=0
    out=1;
else
    out=0;
end
end