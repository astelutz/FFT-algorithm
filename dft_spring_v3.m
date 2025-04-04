% Andrew Lutz - Spring 2025
%% Initialization
close all
clear,clc
file = "fourier_dataset_2.txt"; % grabs our data file
data = readmatrix(file); % converts our data file into a matrix

if mod(size(data,1),2) == 1
    N = size(data,1)-1; % determines the size of our matrices
    T = data(N+1); % picks the period from the final value in our time column
else
    N = size(data,1);
    T = data(N);
end
t = data(1:N,1); % stores the first column (time) as a column vector
f = data(1:N,2); % stores the second column (f(t)) as a column vector
dt = abs(t(1)-t(2)); % determines our delta t from our time matrix
M = 2; % number of wave modes to sweep for
qmax = floor((N-1)/2);
A(qmax,1) = zeros;
C(qmax,1) = zeros;
D(qmax,1) = zeros;
Mag(qmax,1) = zeros;
inv(N,qmax) = zeros;

%Optional Randomness
%{
for j = 1:N
    %f(j) = f(j) + 0.2*rand;
    %f(j) = f(j) + 2*cos(3*2*pi*j*dt);
end
%}
f = f - mean(f);
%% Fourier Transform
for q = 1:qmax
    for n = 1:N    
        A(q) = (2/N) * f(n) * exp((1i*2*pi*q*(n-1))/N) + A(q);
        C(q) = (2/N) * f(n) * cos((2*pi*q*(n-1))/N) + C(q);
        D(q) = (2/N) * f(n) * sin((2*pi*q*(n-1))/N) + D(q);
        %Mag(q) = sqrt((C(q))^2 + (D(q))^2);
        Mag(q) = abs(A(q));
    end
    for n = 1:N
        inv(n,q) = C(q) * cos(2*pi*q*(n-1)/N) + D(q) * sin(2*pi*q*(n-1)/N);
    end
end
F = sum(inv,2);

%% Results
%
figure(1)
plot((1/T):(1/T):(qmax/T),C)
hold on
plot((1/T):(1/T):(qmax/T),D)
hold off
legend('C Coefficient','D Coefficient')
ylabel('Strength of Coefficient')
xlabel('idk yet')
%
x_ticks=0:0.5:qmax;
x_label{floor((qmax+1)*dt),1}=zeros;
x_label{1} = '0';
for j = 2:2*floor((qmax+1)*dt)
    x_label{j} = sprintf('%d\\pi',j-1);
end
xticklabels(x_label);
%}
figure(2)
plot((1/T):(1/T):(qmax/T),Mag)
ylabel('Power')
xlabel('Frequency (kind of)')
xticks(x_ticks)
xticklabels(x_label);
%
figure(3)
plot(t,f)
hold on
plot(t,F)
legend('f(t)','Inv(F(k))')
ylabel('Amplitude (arbitrary units)')
xlabel('Time (s)')
%}