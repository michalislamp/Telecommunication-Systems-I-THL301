%% A.5
%------------------------2------------------------
Px_F_total_2T = 0;
K = 500;
for i=1:K
    % 2-PAM bits generator
    b = (sign(randn(N, 1)) + 1)/2;
    x = bits_to_2PAM(b);

    % Upsample the 2-PAM bit vector
    upsample_x = (1/Ts)*upsample(x, over);

    % Convolution of X and phi to compute the sum
    X_t = conv(upsample_x,phi)*Ts;

    % Fourier Transform
    X_F = fftshift(fft(X_t, Nf)*Ts);
    % Create the periodogram
    num = power(abs(X_F),2);
    Px_F_total_2T = Px_F_total_2T + (num/T_total);
end
Px_F_estimated_2T = Px_F_total_2T/K;
Px_F_theoritical = 1/T*phi_Energy;

figure;
semilogy(f_axis,Px_F_estimated_2T,f_axis,Px_F_theoritical);
legend("Estimated","Theoritical");
title("P_x(F)");
xlim([-Fs/2 Fs/2]);
legend("Average","Theoritical");