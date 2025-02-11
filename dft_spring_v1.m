% Andrew Lutz - Spring 2025
%% Initialization
close all
clear,clc
file = "fourier_dataset_2.txt"; % grabs our data file
data = readmatrix(file); % converts our data file into a matrix

t = data(:,1); % stores the first column (time) as a column vector
f = data(:,2); % stores the second column (f(t)) as a column vector
sz = size(t,1); % determines the size of our matrices
T = t(sz); % picks the period from the final value in our time column
dt = abs(t(1)-t(2));
n = 2; % number of wave modes to sweep for
C(n,1) = zeros; % establishes a set of matrices for coefficients Cn
D(n,1) = zeros;
A(n,1) = zeros;
G(n,1) = zeros; % establishes a constant for an optional exponential
w(n,1) = zeros;
%f = f - sum(f)/sz; subtracts the average of f from f to emphasize peaks

%% Fourier Transform
C0 = (2/T) * sum(f);
for m = 1:n
    w(m,1) = m * (2*pi)/T;
    for j = 1:sz
        C(m) = C(m) + (2/T) * f(j) * cos(w(m,1)*t(j)) * dt;
        Cn(j,m) = C(m);
        disp(C(m))
        D(m) = D(m) + (2/T) * f(j) * sin(w(m,1)*t(j)) * dt;
        Dn(j,m) = C(m);
        A(m) = sqrt((C(m)^2)+(D(m)^2));
        G(m) = G(m) + (2/T) * f(j) * exp(1i*w(m,1)*t(j)) * dt;
    end
end

%% Results
disp('C coefficients')
disp(C)