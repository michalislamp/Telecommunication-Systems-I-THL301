%% A.4

% Creating the pulse
[phi, t] = srrc_pulse(T, over, A_SRRC, a);

% Upsample the 2 sequencies
XIn_up = (1/Ts)*upsample(XIn, over);
XI_up_time = 0:Ts:N*T-Ts;

XQn_up = (1/Ts)*upsample(XQn, over);
XQ_up_time = XI_up_time;

% Convolutions and time
XIt = Ts*conv(XIn_up, phi);
XQt = Ts*conv(XQn_up, phi);
conv_t = t(1) + XQ_up_time(1):Ts:t(end) + XQ_up_time(end);

% Output waveform plots
figure;
plot(conv_t, XIt);
title("Output waveform XI(t)");
xlabel("t");
ylabel("XI(t)");
xlim([conv_t(1) conv_t(end)]);

figure;
plot(conv_t, XQt);
title("Output waveform XQ(t)");
xlabel("t");
ylabel("XQ(t)");
xlim([conv_t(1) conv_t(end)]);

% Frequency Range
f_axis = linspace(-Fs/2,(Fs/2-Fs/Nf), Nf);

% Periodograms creation
XIF = fftshift(fft(XIt, Nf)*Ts);
num = power(abs(XIF),2);
T_total = conv_t(end) - conv_t(1);
PXIF = num/T_total;

XQF = fftshift(fft(XQt, Nf)*Ts);
num = power(abs(XQF),2);
T_total = conv_t(end) - conv_t(1);
PXQF = num/T_total;

% Periodogram plots
figure;
plot(f_axis, PXIF);
title("Periodogram of XI(t)");
xlabel("f");
ylabel("Px(F)");

figure;
plot(f_axis, PXQF);
title("Periodogram of XQ(t)");
xlabel("f");
ylabel("Px(F)");
