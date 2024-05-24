%% A.3
%-----------------------3------------------------
K_options = [1000; 50];
N_options = [200; 50];
figure;
for i=1:length(K_options)
    Px_F_total = 0;
    K_3 = K_options(i);
    N_3 = N_options(i);
    for j=1:K_3
        % 2-PAM bits generator
        b = (sign(randn(N_3, 1)) + 1)/2;
        x = bits_to_2PAM(b);

        % Upsample the 2-PAM bit vector
        upsample_x = (1/Ts)*upsample(x, over);
        upsample_x_time = 0:Ts:N_3*T-Ts;

        % Convolution of X and phi to compute the sum
        X = conv(upsample_x,phi)*Ts;
        X_time = upsample_x_time(1)+t(1):Ts:upsample_x_time(end)+t(end);

        % Create the periodogram
        X_F = fftshift(fft(X, Nf)*Ts);
        num = power(abs(X_F),2);
        T_total = X_time(end) - X_time(1);
        Px_F_total = Px_F_total + (num/T_total);
    end
    Px_F_estimated_3 = Px_F_total/K_3;
    Px_F_theoritical = 1/T*phi_Energy;

    subplot(1,2,i);
    semilogy(f_axis,Px_F_estimated_3);
    hold on;
    semilogy(f_axis,Px_F_theoritical);
    legend("Estimated","Theoritical");
    ylabel("P_x(F)");
    xlabel("Frequency");
    title(["K = ",K_3," N = ",N_3])
    xlim([-Fs/2 Fs/2]);
    hold off;
end