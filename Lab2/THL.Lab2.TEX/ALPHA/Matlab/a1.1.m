clear all;
close all;
%% Variable decleration
T = 0.001;
over = 10;
Ts = T / over;
Fs = 1 / Ts;
A = 4;
a = 0.5;
K = 500;
N = 100;
Nf = 4096;
%% A.1

% Pulse generator
grid on;
hold on;
[phi, t] = srrc_pulse(T, over, A, a);
plot(t,phi);
xlim([-A*T A*T]);
title("SRRC Pulse");
xlabel("Time");
ylabel("\phi(t)");
hold off;

% Frequency Range
f_axis = linspace(-Fs/2,(Fs/2-Fs/Nf), Nf);
% Fourier Transform 
phi_F = fftshift(fft(phi,Nf)*Ts);
abs_phi_F = abs(phi_F);
phi_Energy = power(abs_phi_F,2);
figure;
semilogy(f_axis, phi_Energy);
title("Energy Spectrum of SRRC");
xlabel("Frequency");
ylabel("Logarithmic");
xlim([-Fs/2 Fs/2]);
