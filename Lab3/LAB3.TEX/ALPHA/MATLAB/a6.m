%% A.6

% Input signal creation 
Xmod_t = XImod_t + XQmod_t;

% Plot of Xmod(t)
figure;
plot(conv_t, Xmod_t);
title("Input waveform X(t)mod");
xlabel("t");
ylabel("X(t)");
xlim([conv_t(1) conv_t(end)]);

% Periodogram creation
Xmod_F = fftshift(fft(Xmod_t, Nf)*Ts);
num = power(abs(Xmod_F),2);
T_total = conv_t(end) - conv_t(1);
PXmod = num/T_total;

% Periodogram plot
figure;
plot(f_axis, PXmod);
title("Periodogram of X(t)mod");
xlabel("f");
ylabel("Px(F)");