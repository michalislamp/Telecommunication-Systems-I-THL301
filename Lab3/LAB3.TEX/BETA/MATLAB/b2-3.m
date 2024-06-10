%% B.2 and B.3
figure;
semilogy(SNRdb, P_theor_symbol);
hold on;
semilogy(SNRdb, P_symbol);
legend ("Theoretical","Experimental")
hold off;
grid on;
title("Symbol Comparisson");

figure;
semilogy(SNRdb, P_theor_bit);
hold on;
semilogy(SNRdb, P_bit);
legend ("Theoretical","Experimental")
hold off;
grid on;
title("Bits Comparisson");