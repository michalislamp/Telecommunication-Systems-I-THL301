%% A.10

% Filter the waveforms
XIt_demodulated = Ts*conv(XImod_noise, phi);
XQt_demodulated = Ts*conv(XQmod_noise, phi);
conv_t2 = t(1) + conv_t(1):Ts:t(end) + conv_t(end);

% Plot the filtered waveforms
figure;
plot(conv_t2, XIt_demodulated);
title("Filtered XI(t)")
xlabel("t");
ylabel("XI(t)");
xlim([conv_t2(1) conv_t2(end)]);

figure;
plot(conv_t2, XQt_demodulated);
title("Filtered XQ(t)")
xlabel("t");
ylabel("XQ(t)");
xlim([conv_t2(1) conv_t2(end)]);

XIF_demodulated = fftshift(fft(XIt_demodulated, Nf)*Ts);
num = power(abs(XIF_demodulated),2);
T_total = conv_t2(end) - conv_t2(1);
PXI_demodulated = num/T_total;

XQF_demodulated = fftshift(fft(XQt_demodulated, Nf)*Ts);
num = power(abs(XQF_demodulated),2);
T_total = conv_t2(end) - conv_t2(1);
PXQ_demodulated = num/T_total;

% Plot the periodogram of filtered waveforms
figure;
plot(f_axis, PXI_demodulated);
title("Periodogram of XI(t)")
xlabel("f");
ylabel("PxI(F)");

figure;
plot(f_axis, PXQ_demodulated);
title("Periodogram of XQ(t)")
xlabel("f");
ylabel("PxQ(F)");