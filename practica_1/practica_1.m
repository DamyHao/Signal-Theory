

% Fem la funcio per qualsevol llargada dels vectors:
function y = iir_fd1(a, b, x, N)
% INPUT: N es longitud del vector.
lenghtA = length(a);
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


%Infinite Impulse Response (IIR) system: there is recursive
%component, h[n] has infinite length, the system may be
%unstable (but not necessarily)