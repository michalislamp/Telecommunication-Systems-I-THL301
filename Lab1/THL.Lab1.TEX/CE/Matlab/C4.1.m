%% C.4.1
phi_flip = flip(phi);
Z_conv = conv(X_conv, phi_flip)*Ts;
Z_time = X_time(1)+t(1):Ts:X_time(end)+t(end);

figure
plot(Z_time, Z_conv)
xlabel("Time")
xlim ([Z_time(1) Z_time(end)]) ;
title("Z(t)");
