%% C.3
[phi, t] = srrc_pulse(T, over, A, a);
X_conv = conv(x_delta, phi)*Ts;
X_time = t_pulse(1)+t(1):Ts:t_pulse(end)+t(end);

figure
plot(X_time, X_conv)
xlabel("TIme")
xlim([X_time(1) X_time(end)]) ;
title("X(t)");