function d = dft(f)
    if mod(size(f,1),2) == 1
        N = size(f,1)-1; % determines the size of our matrices
    else
        N = size(f,1);
    end
    qmax = floor((N-1)/2);
    A(qmax,1) = zeros;
    d(qmax,1) = zeros;
    f = f - mean(f);
    for q = 1:qmax
        for n = 1:N    
            A(q) = (2/N) * f(n) * exp((1i*2*pi*q*(n-1))/N) + A(q);
            d(q) = abs(A(q));
        end
    end
end