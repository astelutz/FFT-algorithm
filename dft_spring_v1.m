% Andrew Lutz - Spring 2025
%% Initialization
close all
clear,clc
file = "fourier_dataset.txt"; % grabs our data file
data = readmatrix(file); % converts our data file into a matrix

t = data(:,1); % stores the first column (time) as a column vector
f = data(:,2); % stores the second column (f(t)) as a column vector
sz = size(t,1); % determines the size of our matrices
T = t(sz); % picks the period from the final value in our time column
n = 2; % number of wave modes to sweep for
C(n,1) = zeros; % establishes a set of matrices for coefficients Cn
D(n,1)= zeros;
w(n,1) = zeros;

%% Fourier Transform
C0 = (2/T) * sum(f);
for i = 1:n
    w(i,1) = i * (2*pi);
    for j = 1:sz
        C(i) = (2/T) * (f(j)*cos(w(i,1)*t(j))) + C(i);
        D(i) = (2/T) * (f(j)*sin(w(i,1)*t(j))) + D(i);
    end
end