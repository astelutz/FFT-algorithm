% Andrew Lutz - Spring 2025
%% Initialization
close all
clear,clc
file = "EXPERIMENT.txt"; % data imported as column matrices
data = readmatrix(file);

d = size(data,1); % stores size of our data matrix

%pow2 = nextpow2(size(data,1)); % calculates next highest power of 2 of the
% size of our data matrix, which must be of size 2^n for fft to work.

% An alternative form of determining the next size of 2^n necessary if the
% pow2 command is unavailable, it utilizes bit shifting
Nbit = size(data,1);
numdiv = 0; % this value represents the number of powers of 2 present in
% the value for the size of the matrix
while Nbit >= 1
    Nbit = bitshift(Nbit,-1); %shifting the bit divides by 2 repeatedly
    numdiv = numdiv + 1; % every repeat adds 1 to the counter
end

f = zeros(2^numdiv,1); %zero matrix with zeros padded out past where our
% data will reach

% Optional randomness
%{
for j = 1:size(f,1)
    f(j) = f(j) + 0.5 + 1*rand;
end
%}

N = size(f,1);
h = hann(d); % Hann window to make data smoother, research windowing info
% to understand
f(1:d) = data(1:d,2); %layers our data onto the zero matrix
avg = mean(f(1:d));
f(1:d) = f(1:d) - avg; % subtracts mean value of f from matrix
f(1:d) = f(1:d) .* h; % applies our Hann window
dt = abs(data(1,1)-data(2,1)); % dt necessary for axes plotting later on in
% the code, check your time columns before doing anything

%% Outputs
x = mfft(f);
magx = abs(x); % x is a complex number, so the magnitude of x is the
% quantity we actually care about
y = normalize(magx,'range');
z = y ./ (avg/100);
%% Plotting
plot((1/(dt*N)*(0:N/2-1)),z(1:N/2))
x_ticks=0:dt:N/2;
x_label{floor((d+1)*dt),1}=zeros;
x_label{1} = '0';
for j = 2:2*floor((d+1)*dt) % this loop will create our n*pi x axis
    x_label{j} = sprintf('%d\\pi',j-1);
end
xticklabels(x_label);
%}
ylabel('Power')
xlabel('Frequency (kind of)')
%% FFT Function
function V = mfft(f) % declares the FFt function
    N = size(f,1);
    W = exp(-2*pi*1i/N); % establishes roots of unity based on matrix size,
    % this is how the code "sweeps the frequencies"
    if N == 1 
        V = f;
        return
    else
        fe(1:N/2,1) = zeros; % establishes matrix of even indexed f values
        fo(1:N/2,1) = zeros; % establishes matrix of odd indexed f values
        for n = 1:(N/2)
            fe(n) = f(2*n-1);   
            fo(n) = f(2*n);
        end
        % recursive part of FFT below, we further divide & conquer the
        % until it reaches size 1
        Ve = mfft(fe); 
        Vo = mfft(fo);
%-------------------------------------------------------
        V(1:N,1) = zeros;
        for i = 1:(N/2)
            V(i) = Ve(i) + (W^(i-1)) * Vo(i); % applying roots of unity
            V(N/2 + i) = Ve(i) - (W^(i-1)) * Vo(i);
        end
    end
end