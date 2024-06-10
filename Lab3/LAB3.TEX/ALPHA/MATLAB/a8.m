%% A.8

% Noise creation
SNR = 10;
SNR = 20;
variance = (10*(A^2))/(Ts*(10^(SNR/10)));
gaussian_noise = sqrt(variance)*randn(1,length(Xmod_t));

% Adding noise to the waveform