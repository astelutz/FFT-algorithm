function d = idft(Re, Im)
    N = 2 * (size(Re,1) + 1);
    qmax = floor((N-1)/2);
    inv(N,qmax) = zeros;
    for q = 1:qmax
        for n = 1:N
            inv(n,q) = Re(q) * cos(2*pi*q*(n-1)/N) + Im(q) * sin(2*pi*q*(n-1)/N);
        end
    end
    d = sum(inv,2);
end