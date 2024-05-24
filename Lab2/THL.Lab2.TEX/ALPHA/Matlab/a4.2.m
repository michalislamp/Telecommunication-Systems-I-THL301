%% A.4
%---------------------2----------------------
figure;
subplot(1,2,1);
plot(f_axis, Px_F_estimated);
hold on;
plot(f_axis,Px_F_estimated_4PAM)
title("P_X(F)");
xlabel("Frequency");
legend("2PAM","4PAM");
xlim([-Fs/2 Fs/2]);
hold off;

subplot(1,2,2);
semilogy(f_axis,Px_F_estimated);
hold on;
semilogy(f_axis,Px_F_estimated_4PAM);
title("P_X(F)");
xlabel("Frequency");
legend("2PAM","4PAM");
xlim([-Fs/2 Fs/2]);
hold off;