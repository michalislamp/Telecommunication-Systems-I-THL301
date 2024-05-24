%% A.4
%----------------------1-----------------------
P_x_4PAM_total = 0;
for i=1:K
    % 2-PAM bits generator
    b = (sign(randn(N, 1)) + 1)/2;
    x = bits_to_4PAM(b);

    % Upsample the 2-PAM bit vector
    x_up = (1/Ts)*upsample(x, over);
    x_up_time = 0:Ts:(N/2)*T-Ts;

    % Convolution of X and phi to compute the sum
    X = conv(x_up,phi)*Ts;
    X_time = x_up_time(1)+t(1):Ts:x_up_time(end)+t(end);

    % Fourier transform of X(t)
    x_F_4PAM = fftshift(fft(X,Nf)*Ts);

    % Create the periodogram
    T_total = X_time(end) - X_time(1);
    P_x_4PAM_num = power(abs(x_F_4PAM),2);
    P_x_4PAM_total = P_x_4PAM_total + (P_x_4PAM_num/T_total);
end

Px_F_estimated_4PAM = P_x_4PAM_total/K;
Px_F_theoritical_4PAM = 5/T*phi_Energy;

figure;
semilogy(f_axis,Px_F_estimated);
hold on;
semilogy(f_axis,Px_F_theoritical);
legend("Estimated","Theoritical");
title("P_x(F)");
xlim([-Fs/2 Fs/2]);
hold off;