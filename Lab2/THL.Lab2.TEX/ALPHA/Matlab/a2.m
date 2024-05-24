%% A.2

% 2-PAM bits generator
b = (sign(randn(N, 1)) + 1)/2;
x = bits_to_2PAM(b);

% Upsample the 2-PAM bit vector
upsample_x = (1/Ts)*upsample(x, over);
upsample_x_time = 0:Ts:N*T-Ts;

% Convolution of X and phi to compute the sum
X = conv(upsample_x,phi)*Ts;
X_time = upsample_x_time(1)+t(1):Ts:upsample_x_time(end)+t(end);
figure;
plot(X_time, X);
xlabel("Time");
ylabel("X(t)");
xlim([X_time(1) X_time(end)]);