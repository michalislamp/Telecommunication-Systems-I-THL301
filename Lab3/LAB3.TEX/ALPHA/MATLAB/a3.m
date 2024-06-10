%% A.3

% For bits in the [0,2N]
bitXI = bits(1:2*N);
XIn = bits_to_4_PAM(bitXI, A);

% For bits in the [2N+1,4N]
bitXQ = bits(2*N+1:4*N);
XQn = bits_to_4_PAM(bitXQ, A);