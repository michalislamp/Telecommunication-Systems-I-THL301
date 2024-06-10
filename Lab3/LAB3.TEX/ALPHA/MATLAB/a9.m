%% A.9

% Demodulating the signal
XImod_noise = Xmod_noise.*cos(2*pi*Fo*conv_t);
XQmod_noise = -Xmod_noise.*sin(2*pi*Fo*conv_t);

% Plot the signals
figure;
plot(conv_t, XImod_noise);
title("XI(t)cos(2piFot)")
xlabel("t");
ylabel("XI(t)");
xlim([conv_t(1) conv_t(end)]);

figure;
plot(conv_t, XQmod_noise);
title("-XQ(t)sin(2piFot)")
xlabel("t");
ylabel("XQ(t)");
xlim([conv_t(1) conv_t(end)]);

% Periodograms creation
XImod_noise_F = fftshift(fft(XImod_noise, Nf)*Ts);
num = power(abs(XImod_noise_F),2);
T_total = conv_t(end) - conv_t(1);
PXImod_noise = num/T_total;

XQmod_noise_F = fftshift(fft(XQmod_noise, Nf)*Ts);
num = power(abs(XQmod_noise_F),2);
T_total = conv_t(end) - conv_t(1);
PXQmod_noise = num/T_total;

% Plot the periodogram
figure;
plot(f_axis, PXImod_noise);
title("Periodogram of XI(t)")
xlabel("f");
ylabel("Px(F)");

figure;
plot(f_axis, PXQmod_noise);
title("Periodogram of XQ(t)")
xlabel("f");
ylabel("Px(F)");