%% A.15

bitXI_received = PAM_4_to_bits(detected_XI, A);
bitXQ_received = PAM_4_to_bits(detected_XQ, A);

estimated_bits = [bitXI_received bitXQ_received];
difference = estimated_bits - bits;
errors_in_bits = 0;
for i = 1:4*N
        if (difference(i) ~= 0)
            errors_in_bits = errors_in_bits + 1;
        end
end
errors_in_bits