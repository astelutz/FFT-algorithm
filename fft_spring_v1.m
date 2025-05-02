% Andrew Lutz - Spring 2025
%% Initialization
close all
clear,clc
%file = "fourier_dataset_2.txt";
file = "experiment_5_HG.txt";
%file = "sample-data.txt";
data = readmatrix(file);

%pow2 = nextpow2(size(data,1));

d = size(data,1);
Nbit = size(data,1);
numdiv = 0;
while Nbit >= 1
    Nbit = bitshift(Nbit,-1);
    numdiv = numdiv + 1;
end

f = zeros(2^numdiv,1);
%{
for j = 1:size(f,1)
    f(j) = f(j) + 0.5 + 1*rand;
end
%}
N = size(f,1);
h = hann(d);
f(1:d) = data(1:d,2);
avg = mean(f(1:d));
f(1:d) = f(1:d) - mean(f(1:d));
f(1:d) = f(1:d) .* h;
dt = abs(data(1,1)-data(2,1));

%% Outputs
%W = exp(2*pi*1i/N);
%Nt = N;
x = mfft(f);
magx = abs(x);
y = normalize(magx,'range');
z = y ./ (avg/100);
%% Plotting
plot((1/(dt*N)*(0:N/2-1)),z(1:N/2))
%
x_ticks=0:dt:N/2;
x_label{floor((d+1)*dt),1}=zeros;
x_label{1} = '0';
for j = 2:2*floor((d+1)*dt)
    x_label{j} = sprintf('%d\\pi',j-1);
end
xticklabels(x_label);
%}
ylabel('Power')
xlabel('Frequency (kind of)')
%% FFT Function
function V = mfft(f)
    N = size(f,1);
    W = exp(-2*pi*1i/N);
    if N == 1
        V = f;
        return
    else
        fe(1:N/2,1) = zeros;
        fo(1:N/2,1) = zeros;
        for n = 1:(N/2)
            fe(n) = f(2*n-1);
            fo(n) = f(2*n);
        end
        Ve = mfft(fe);
        Vo = mfft(fo);
%-------------------------------------------------------
        V(1:N,1) = zeros;
        for i = 1:(N/2)
            V(i) = Ve(i) + (W^(i-1)) * Vo(i);
            V(N/2 + i) = Ve(i) - (W^(i-1)) * Vo(i);
        end
    end
end
%%%%