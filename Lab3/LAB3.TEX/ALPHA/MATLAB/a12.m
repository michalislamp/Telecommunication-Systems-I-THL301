%% A.12

% Estimating the symbols value
for i = 1:N
    detected_XI(i) = detect_4_PAM(XIt_demodulated_sampling(i), A);
    detected_XQ(i) = detect_4_PAM(XQt_demodulated_sampling(i), A);
end