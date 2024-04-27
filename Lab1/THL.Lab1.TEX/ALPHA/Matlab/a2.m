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