% Andrew Lutz - Spring 2025
%% Initialization
close all
clear,clc
file = "fourier_dataset_2.txt"; % grabs our data file
data = readmatrix(file); % converts our data file into a matrix

t = data(:,1); % stores the first column (time) as a column vector
f = data(:,2); % stores the second column (f(t)) as a column vector
N = size(t,1); % determines the size of our matrices
T = t(N); % picks the period from the final value in our time column
dt = abs(t(1)-t(2)); % determines our delta t from our time matrix
M = 2; % number of wave modes to sweep for
qmax = floor((N/2)-1);
A(qmax,M) = zeros;
B(qmax,M) = zeros;
C(qmax,M) = zeros;

%% Fourier Transform
for q = 1:qmax
    for n = 1:N
        A(q) = (2/N) * f(n) * exp((1i*2*pi*q*n)/N) + A(q);
        B(q) = (2/N) * f(n) * cos((2*pi*q*n)/N) + B(q);
        C(q) = (2/N) * f(n) * sin((2*pi*q*n)/N) + C(q);
    end
    disp('C Coeffcient')
    disp(B(q))
    disp('D Coefficient')
    disp(C(q))
end

%% Results
figure(1)
plot(t,B)