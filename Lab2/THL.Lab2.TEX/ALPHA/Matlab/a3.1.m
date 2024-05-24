%% A.3
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
