%% A.5
T = 2*T;
over = 2*over;
% Create the new pulse
[phi, t] = srrc_pulse(T, over, A, a);

phi_F = fftshift(fft(phi, Nf)*Ts);
abs_phi_F = abs(phi_F);
phi_Energy = power(abs_phi_F,2);

% 2-PAM bits generator
b = (sign(randn(N, 1)) + 1)/2;
x = bits_to_2PAM(b);

% Upsample the 2-PAM bit vector
upsample_x = (1/Ts)*upsample(x, over);
upsample_x_time = 0:Ts:N*T-Ts;

% Convolution of X and phi to compute the sum
X = conv(upsample_x,phi)*Ts;
X_time = upsample_x_time(1)+t(1):Ts:upsample_x_time(end)+t(end);

%-----------------------1-----------------------
% Create the periodogram
X_F = fftshift(fft(X, Nf)*Ts);
num = power(abs(X_F),2);
T_total = X_time(end) - X_time(1);
Px_F = num/T_total;

% Visualize Px(F) and RF Rx(F) using log scale
figure;
subplot(1,2,1);
plot(f_axis, Px_F);
title("Px(F)");
xlabel("Frequency");
xlim([-Fs/2 Fs/2]);

subplot(1,2,2);
semilogy(f_axis,Px_F);
title("Px(F)");
xlabel("Frequency");
ylabel("Logarithmic");
xlim([-Fs/2 Fs/2]);
