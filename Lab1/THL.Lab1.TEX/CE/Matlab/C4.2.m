%% C.4.2
% Sampling Z(t) 
index_first = find(abs(Z_time) < T);
index_last = find(abs(Z_time - t_pulse(end)) < T);

Z_sampled = Z_conv(index_first : index_last);
Z_sampled = downsample(Z_sampled, over);

figure
hold on;
plot(Z_time, Z_conv);
stem(0:T:(N-1)*T, Z_sampled, "k");
legend ("Z(t)", "Z_s(t)");
xlim([Z_time(1) Z_time(end)]);
