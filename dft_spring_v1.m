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
A(n,1) = zeros; % establishes a constant for an optional exponential
w(n,1) = zeros;
g(T,1) = zeros;
%f = f - sum(f)/sz; %subtracts the average of f from f to emphasize peaks

%% Fourier Transform
C0 = (2/T) * sum(f);
for m = 1:n
    w(m,1) = m * (2*pi)/T; % omega
    for j = 1:T
        C(m) = C(m) + (2/T) * f(j) * cos((w(m)*(j-1))); % C coefficient
        D(m) = D(m) + (2/T) * f(j) * sin(w(m)*(j-1)); % D coefficient
        A(m) = A(m) + (2/T) * f(j) * exp(1i*w(m)*(j-1)); % exponential
        g(j,1) = C0 + C(m) * cos(w(m)*(j-1)*T) + D(m) * sin(w(m)*(j-1)*T);
    end
end


%% Results
disp('C coefficients')
disp(C)
disp('D coefficients')
disp(D)

figure(1)
plot(g)

figure(2)
plot(f)