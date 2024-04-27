%% A.3.

c1 = zeros(1,length(f_axis))+T/10^3;  %Limit C1 = T/10^3
c2 = zeros(1,length(f_axis))+T/10^5;  %Limit C2 = T/10^5

figure(4);
for i=1:length(a)
    semilogy(f_axis, phi_energy{i});
    hold on;
end
grid on;
title('Energy Spectrums of SRRC');
xlabel('Frequency');
ylabel('Logarithmic');
plot(f_axis, c1);
plot(f_axis, c2);
legend("a=0", "a=0.5","a=1", "c1 = T/10^3","c2 = T/10^5");
xlim([-Fs/2 Fs/2])
hold off;
