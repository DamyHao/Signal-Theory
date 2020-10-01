%Infinite Impulse Response (IIR) system: there is recursive
%component, h[n] has infinite length, the system may be
%unstable (but not necessarily)

clear all;
clc;
close all;

a=[1 -7/8 3/32];
b=[1 -2/3];

N = 10;
for k = 0:N-1
    x(k+1)=delta(k)-2*delta(k-1)+4*delta(k-2);
end

 y = iir_fd1(a, b, x, N);
 z = iir_fd1M(a, b, x, N);
 z2 = iir_fd2M(a, b, x, N);
 
 

 figure();
 %Per valors discrets
stem(1:10,x,'-o')
grid on

figure();
stem(1:10,z,'-o');
figure();
 stem(1:10,z2,'-o');
 

function out = delta(in)
out = 0;
if(in == 0)
    out = 1;
end
end

% Fem la funcio per qualsevol llargada dels vectors:
function y = iir_fd1(a, b, x, N)
% INPUT: N es longitud del vector.
lengthA = length(a);
lengthB = length(b);
sum = 0;

y = zeros( 1, N);

% calculant en punt n
for n = 1:N
    for it = 0:(lengthA-1)
        if (n-it) > 0
            sum = sum + x(n-it)*a(it+1);
        end
    end
    
    sum2 = 0;
    %Comença al 1 perq no tenim acces a
    for itt = 1:lengthB
        if (n-it) > 1
            sum2 = sum2 + y(n-it)*b(it);
        end
    end
    y(n) = sum - sum2;
end

end

function y=iir_fd1M(a,b,x,N)
y(1)=b(1)*x(1);
y(2)=b(1)*x(2)+b(2)*x(1)-a(2)*y(1);
for k=1:N-2
    compox=b(1)*x(k+2)+b(2)*x(k+1);
    compoy=-a(2)*y(k+1)-a(3)*y(k);
    y(k+2)=(1/(a(1)))*(compox+compoy);
end
end

function y=iir_fd2M(a,b,x,N)
y(1)=b(1)*x(1);
y(2)=b(1)*x(2)+b(2)*x(1)-a(2)*y(1);
for k=1:N-2
    y(k+2)=(1/(a(1)))*(b(1)*x(k+2)+b(2)*x(k+1)-a(2)*y(k+1)-a(3)*y(k));
end
end
