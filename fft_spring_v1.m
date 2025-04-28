% Andrew Lutz - Spring 2025
%% Initialization
close all
clear,clc
file = "fourier_dataset_2.txt";
data = readmatrix(file);

%pow2 = nextpow2(size(data,1));

Nbit = size(data,1);
numdiv = 0;
while Nbit >= 1
    Nbit = bitshift(Nbit,-1);
    numdiv = numdiv + 1;
end

f = zeros(2^numdiv,1);
N = size(f,1);
f(1:size(data,1)) = data(1:size(data,1),2);
%% D&C and Transform
%W = exp(2*pi*1i/N);
%Nt = N;
x = mfft(f);
magx = abs(x);
plot(magx)
function V = mfft(f)
    N = size(f,1);
    Nt = N;
    W = exp(2*pi*1i/N);
    if Nt == 1
        V = f;
        return
    else
        fe(1:Nt/2,1) = zeros;
        fo(1:Nt/2,1) = zeros;
        for n = 1:(Nt/2)
            fe(n) = f(2*n-1);
            fo(n) = f(2*n);
        end
        %Nt = Nt/2;
        Ve = mfft(fe);
        Vo = mfft(fo);
%-------------------------------------------------------
        V(1:Nt,1) = zeros;
        for i = 1:(Nt/2)
            V(i) = Ve(i) + (W^(i-1)) * Vo(i);
            V(Nt/2 + i) = Ve(i) - (W^(i-1)) * Vo(i);
        end
    end
end
%%%%