%% A.5

Fo = 200;
XImod_t = 2*XIt.*cos(2*pi*Fo*conv_t);
XQmod_t = -2*XQt.*sin(2*pi*Fo*conv_t);

% Plots of the new waveforms
figure;
plot(conv_t, XImod_t);
title("2XI(t)cos(2piFot)");
xlabel("t");
ylabel("XI(t)");
xlim([conv_t(1) conv_t(end)]);

figure;
plot(conv_t, XQmod_t);
title("-2XQ(t)sin(2piFot)");
xlabel("t");
ylabel("XQ(t)");
xlim([conv_t(1) conv_t(end)]);

% Periodograms creation
XIFmod = fftshift(fft(XImod_t, Nf)*Ts);
num = power(abs(XIFmod),2);
T_total = conv_t(end) - conv_t(1);
PXIFmod = num/T_total;

XQFmod = fftshift(fft(XQmod_t, Nf)*Ts);
num = power(abs(XQFmod),2);
T_total = conv_t(end) - conv_t(1);
PXQFmod = num/T_total;

% Periodogram plots
figure;
plot(f_axis, PXIFmod);
title("Periodogram of XI(t)mod");
xlabel("f");
ylabel("Px(F)");

figure;
plot(f_axis, PXQFmod);
title("Periodogram of XQ(t)mod");
xlabel("f");
ylabel("Px(F)");