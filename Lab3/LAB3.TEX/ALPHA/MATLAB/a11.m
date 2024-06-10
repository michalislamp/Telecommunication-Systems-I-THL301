%% A.11

XIt_demodulated_sampling= XIt_demodulated((2*A_SRRC*T/Ts)+1:over:
    length(XIt_demodulated)-(2*A_SRRC*T/Ts));

XQt_demodulated_sampling= XQt_demodulated((2*A_SRRC*T/Ts)+1:over:
    length(XQt_demodulated)-(2*A_SRRC*T/Ts));

for i=1:N
    sampling(i,1)=XIt_demodulated_sampling(i); 
    sampling(i,2)=XQt_demodulated_sampling(i);
end
scatterplot(sampling)