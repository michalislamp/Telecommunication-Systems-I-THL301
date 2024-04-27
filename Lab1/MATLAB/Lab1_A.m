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

%% A.2.

Fs = 1/Ts;  %Sampling Period
%Nf = 1024;
Nf = 2048;  %Step
phi_F = {}; %Initialization of Fourier
phi_energy = {}; %Initialization

f_axis =  linspace(-Fs/2,(Fs/2-Fs/Nf),Nf);  %Freq range
figure(2)
hold on;
grid on;

%Set energy spectrum and plot usinf fft,fftshift
for i=1:length(a)
        phi_F{i} = fftshift(fft(phi{i}, Nf)*Ts);
        phi_energy{i} = power(abs(phi_F{i}),2);
        plot(f_axis, phi_energy{i});
end
legend("a=0", "a=0.5","a=1");
title('Energy Spectrums of SRRC(Nf=2048)');
xlabel('Frequency');
xlim([-Fs/2 Fs/2])
hold off;

%Set energy spectrum using semilogy
figure(3)
for i=1:length(a)
    semilogy(f_axis, phi_energy{i});
    hold on;
end
grid on;
legend("a=0", "a=0.5","a=1");
xlim([-Fs/2 Fs/2])
title('Energy Spectrums of SRRC(Nf=2048)');
xlabel('Frequency');
ylabel('Logarithmic');
hold off;

%% A.3.

c1 = zeros(1,length(f_axis))+T/10^3;  %Limit C1 = T/10^3
c2 = zeros(1,length(f_axis))+T/10^5;  %Limit C2 = T/10^5

figure(4);
for i=1:length(a)
    semilogy(f_axis, phi_energy{i});
    hold on;
end
grid on;
title('Energy Spectrums of SRRC');
xlabel('Frequency');
ylabel('Logarithmic');
plot(f_axis, c1);
plot(f_axis, c2);
legend("a=0", "a=0.5","a=1", "c1 = T/10^3","c2 = T/10^5");
xlim([-Fs/2 Fs/2])
hold off;











