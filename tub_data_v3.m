% Andrew Lutz - 6/16/2022
close all
clear,clc
% Below are the variables to change in this analysis.

%file_num = 12;   % The number of the experiment file.
file = "Exp57tilt.csv";
fWin = 150;     % The fourier time window.
Win_min = 1;     % The time window to begin analysis at.
Win_max = 901;
mode = 3;
tWin = fWin*0.1;
sigma = pi/3;               % Value for our sigma.
k = fWin/2-1;               % Maximum number of q values usable for fWin.
Nsensor = 1;
L = 0.74; %length of the tub
depth = 0.15; %depth of the tub
y = 0.02; %tilt height
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
switch file_num
    case 5
        file = "Experiment5.xlsx";
        subtitle1 = "Pool Data 2020 - Experiment 5";
    case 8
        file = "Experiment8.xlsx";
        subtitle1 = "Pool Data 2020 - Experiment 8";
    case 11
        file = "Experiment11.xlsx";
        subtitle1 = "Pool Data 2020 - Experiment 11";
    case 12
        file = "Experiment12.xlsx";
        subtitle1 = "Pool Data 2020 - Experiment 12";
end % Chooses one of the three experiment files to analyze.
%}
subtitle1 = "Deep Tub Data 2023 - (Depth = 15cm, Y = 10cm)";

switch tWin
    case 15
        subtitle2 = "(Fourier Window = 15 Seconds)";
    case 30
        subtitle2 = "(Fourier Window = 30 Seconds)";
    case 45
        subtitle2 = "(Fourier Window = 45 Seconds)";
    case 60
        subtitle2 = "(Fourier Window = 60 Seconds)";
    case 90
        subtitle2 = "(Fourier Window = 90 Seconds)";
    case 120
        subtitle2 = "(Fourier Window = 120 Seconds)";
end % Chooses a time window in which to analyze
    % and sets graph subtitles for future use.
EXPE = readmatrix(file);       % Accesses the values from the experiment file.
T = EXPE(:,1);              % Creates an array of time values.
Tfix = T(Win_min:Win_max-fWin,:);  % Chooses a set amount of those time values.
for i = Win_min:Win_max-fWin
S = EXPE(i:fWin+i,Nsensor+1);       % Stores the values our sensor detected.
SAVG = (sum(S))/fWin;       % Finds the average sensor value.
S = S - SAVG;               % Stores S as new values.
    for q = 1:k;
    B(q) = 0;
    A(q) = 0;
        for n = 1:fWin
        B(q) = (2/fWin)*S(n)*cos((2*pi*q*n)/fWin) + B(q);
        A(q) = (2/fWin)*S(n)*sin((2*pi*q*n)/fWin) + A(q);
        C(q) = sqrt((A(q).^2)+(B(q).^2));
        % These formulas run our data through a fourier analysis.
        end
    end
    D(i,:) = C;             % Stores our fourier values in an array.
end
Dfix = D(Win_min:Win_max-fWin,:);  % Gives us a set range of sensor values
                            % based on 2 of our changing variables.
f = (1:(fWin/2-1))/tWin;    % Calculates the frequencies of our data.

imagesc(f,Tfix,Dfix)
hcb = colorbar
hcb.Title.String = "Amplitude (cm)"
xlabel('Frequency (Hz)', 'FontSize', 10)
xlim([0,max(f)])
ylabel('Time (sec)','FontSize',10)
title('Time vs Frequency')
subtitle([subtitle1 ' ' subtitle2])
% Graphs Time vs Frequency according to our sensor values and time window.
figure(2)
N = mean(Dfix);
M = size(N);
plot(1:M(2),N(1,:));
title('Average Amplitude vs Column Number')
subtitle([subtitle1 '' subtitle2])
ylabel('Average Amplitude (cm)')
xlabel('Column Number')


%Determines best fit frequency colomn from D_fix and calculates the rate of
% decay and plots an approximated exponential decay rate.
Dfix_avg = mean(Dfix);     %Sets a an array of the average of all the amplitudes
qavg = Dfix_avg;           %Sets a dummmy varible to remove data from
data_size=size(Dfix);            %Stores the size of data being used

quad_set = 6;


for h=1:mode   %Calculates number of mode depending on user input
    col_test=false;         %Sets column varible to confirm bad/good column
    while col_test==false   %Checks each column on whether it is increasing or decreasing
        q(:,h)=max(qavg);                             %Finds maxium amplitude value
        [row(h), col_best(h)]=find(q(:,h)==qavg);     %Finds what column maxim amplitude value is located in
        for i=1:quad_set    %Grabs quadrant_set data and checks deterivative sum and sees if it is above 0 or not.
            quad(i,1)=sum(diff(Dfix(i:i+(data_size(1)/quad_set),col_best(h))));    %Calculates deterivative sum
            i=i+(data_size(1)/quad_set)-1;                                                    %Increase quadrant data set capture
        end
        if (sum(quad(:,1)))<0     %Checks if column values are increasing or decreasing
            col_test=true;        %Changes col_test varible to true to break while loop if column is decreasing
       
        else
            if col_best(h)==1     %Checks if best fit column is column 1 *Only needed for experiment 11*
               col_best(h)=2;     %Assigns best fit column to 2 if column 1
               
            end
            qavg(:,col_best(h)-1:col_best(h)+1)=0;   %Deletes bad columnn and surrounding data to eliminate noise
        end
    end
    col_range=round(fWin/60);               %Sets column range
    col_min=col_best(h)-col_range;    %Deletes surrounding data to allow next mode to be located if any...
    col_max=col_best(h)+col_range;    %
    if col_min == 0
        col_min = 1;
    end
                                      %
    qavg(:,col_min:col_max)=0;        %
   
    act_mode(:,h)=Dfix(:,col_best(h));               %Grabs the best fit amplitude and stores it in varible F
    G(:,h)=log(act_mode(:,h)./act_mode(1,h));            %Calculates the linear function
    Gfit(:,h)=polyfit(Tfix,G(:,h),1);           %Uses polyfit to approximate the slope of the linear funciton and the y-intercept
    alpha(:,h)=1/Gfit(1,h);                               %Stores the slope as alpha which is rate of decay of function F
    H(:,h)=act_mode(1,h)*exp((Tfix-min(Tfix))/alpha(:,h));  %Approximates the exponential rate of decay and stores it in function H
    fmode(h) = f(:,col_best(h));
end

%{
figure(3)
plot(Tfix, act_mode)
xlabel('Time (sec)')
ylabel('Amplitude (cm)')
title('Amplitude vs Time (Actual)')
subtitle([subtitle1 ' ' subtitle2])
% Plots the actual values of our amplitude vs our time.

figure(4)
plot(Tfix, H)
xlabel('Time (sec)')
ylabel('Amplitude (cm)')
title('Amplitude vs Time (Approximated)')
subtitle([subtitle1 '' subtitle2])

figure(5)
plot(Tfix, act_mode, 'bo', Tfix, H, 'r')
f5=tiledlayout(3,1,'TileSpacing','Compact');
title('Amplitude vs Time (Approximated and Actual)')
subtitle([subtitle1 ' ' subtitle2])
ylabel('Amplitude (cm)')
xlabel('Time (sec)')
 for c=1:mode   %Determine how many graphs to make depending on mode
     if c==3
         colord='go';    %Sets mode 1 actual to blue dots
         colorl='m';     %Sets mode 1 approximated to magenta line
         m=1;
     elseif c==2
         colord='ro';    %Sets mode 4 actual to red dots
         colorl='y';     %Sets mode 4 approximated to yellow line
         m=4;
     else
         colord='bo';    %Sets mode 3 actual to green dots
         colorl='k';     %Sets mode 3 approximated to black line
         m=3;
     end
     nexttile
     plot(Tfix,act_mode(:,c),colord,Tfix,H(:,c),colorl)
     title(['Mode ',num2str(m)])
 end

figure(6)
plot(Tfix, G)
ylabel('ln(F/F0)')
xlabel('Time (sec)')
title('ln(F/F0) vs Time')
subtitle([subtitle1 ' ' subtitle2])
% Our decay graph
%}
figure(7)
plot(EXPE(:,1),EXPE(:,2))
ylabel('Amplitude')
xlabel('Time (sec)')
title('Amplitude vs Time')
subtitle([subtitle1 ' ' subtitle2])


grav = 9.81;
freq = sort(fmode);

SW = zeros(1,mode);
GW = zeros(1,mode);

for z = 1:mode
    lam = (2*L)/z;
    Wn = (2*pi)/lam;
    SW(1,z) = (z/(2*L))*sqrt(grav*depth);
    GW(1,z) = sqrt(Wn*grav*tanh(Wn*depth))/(2*pi);
    actw(1,z) = freq(1,z);
end

figure(8)
plot(SW,'b')
hold on
plot(GW,'r')
hold on
plot(GW,'ro')
hold on
plot(SW,'bo')
hold on
plot(actw,'mo')
hold on
plot(actw,'m')
legend('Shallow Wave','General Wave','','','Data')
xlabel('Mode')
ylabel('Frequency')

figure(9)
semilogy(f,((N(1,:).*N(1,:)*0.5)/f(1)));
xlabel('Frequency (Hz)')
ylabel('Energy Density (cm^2 / Hz)')
