%% Variable decleration
N = 200;
A = 1;
A_SRRC = 4;
a = 0.5;
T = 0.01;
over = 10;
Ts = T/over;
Fs = 1/Ts;
Nf = 2048;

%% A.1
bits = (sign(randn(4*N, 1)) + 1)/2;