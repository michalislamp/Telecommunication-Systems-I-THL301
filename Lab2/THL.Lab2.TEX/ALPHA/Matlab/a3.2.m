%% A.3
%------------------------2------------------------
Px_F_total = 0;
for i=1:K
    % 2-PAM bits generator
    b = (sign(randn(N, 1)) + 1)/2;
    x = bits_to_2PAM(b);

    % Upsample the 2-PAM bit vector
    upsample_x = (1/Ts)*upsample(x, over);
    %upsample_x_time = 0:Ts:N*T-Ts;

    % Convolution of X and phi to compute the sum
    X = conv(upsample_x,phi)*Ts;
    %X_time = upsample_x_time(1)+t(1):Ts:upsample_x_time(end)+t(end);

    % Create the periodogram
    X_F = fftshift(fft(X, Nf)*Ts);
    num = power(abs(X_F),2);
    T_total = X_time(end) - X_time(1);
    Px_F_total = Px_F_total + (num/T_total);
end
Px_F_estimated = Px_F_total/K;
Px_F_theoritical = 1/T*phi_Energy;

figure;
semilogy(f_axis,Px_F_estimated);
hold on;
semilogy(f_axis,Px_F_theoritical);
legend("Estimated","Theoritical");
title("P_x(F)");
xlim([-Fs/2 Fs/2]);
hold off;