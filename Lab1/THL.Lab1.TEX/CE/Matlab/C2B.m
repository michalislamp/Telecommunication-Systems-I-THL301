%% C.2.B
x_delta = (1/Ts)*upsample(x, over);
t_pulse = 0:Ts:T*N-Ts;

figure
hold on;
stem(t_pulse, x_delta) ;
title("Upscaled Pulse Train")
xlim([0 T*N]) ;
xlabel("Time") ;
