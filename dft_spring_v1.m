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
dt = abs(t(1)-t(2)); % determines our delta t from our time matrix
n = 2; % number of wave modes to sweep for
C(n,1) = zeros; % establishes a set of matrices for coefficients Cn
D(n,1) = zeros; % establishes a set of matrices for coefficients Dn
A(n,1) = zeros; % establishes a constant for an optional exponential
w(n,1) = zeros; % establishes a matrix for our angular frequencies
g(sz-1,n) = zeros; % 
%f = f - sum(f)/sz; %subtracts the average of f from f to emphasize peaks
%% Fourier Transform
C0 = (2/T) * sum(f);
for m = 1:n
    w(m,1) = m * (2*pi); % omega
    for j = 1:sz-1
        C(m) = C(m) + (2/(sz-1)) * f(j) * cos((w(m)/T)*(j-1)); % C coefficient
        D(m) = D(m) + (2/(sz-1)) * f(j) * sin((w(m)/T)*(j-1)); % D coefficient
        A(m) = A(m) + (2/(sz-1)) * f(j) * exp(1i*(w(m)/T)*(j-1)); % exponential
    end
    for j = 1:sz-1
        g(j,m) = C(m) * cos(w(m)*t(j)) + D(m) * sin(w(m)*t(j));
    end
end 
s = sum(g,2);% + C0;
s_sz = size(s,1);

%% Results
disp('C coefficients')
disp(C)
disp('D coefficients')
disp(D)

figure(1)
plot(t(1:s_sz),f(1:s_sz),'r.')
hold on
plot(t(1:s_sz),s,'b')
legend('f(t)','Inv f(w)')
