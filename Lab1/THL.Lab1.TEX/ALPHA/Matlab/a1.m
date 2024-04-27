clear all;
close all;

%% Variable Decleration
T = 0.01;                 %Period
over = 10;                %Oversampling factor
Ts = T/over;              %Sampling Period
A = 4 ;                   %Half duration of the pulse in symbol periods
a = [0, 0.5, 1];          %Roll-off Factor
phi = {};                 %Initialization of truncated SRRC pulse
t = 0;                    %Initialization of time
%% A.1.

figure(1)
grid on;
hold on;

%Set variables for each 'a' and plot
for i=1:length(a)
    [phi{i}, t] = srrc_pulse(T, over, A, a(i));
    plot(t,phi{i});
end
legend("a=0", "a=0.5","a=1")
xlim([-A*T A*T])
title("SRRC Pulses")
xlabel("Time")
hold off;