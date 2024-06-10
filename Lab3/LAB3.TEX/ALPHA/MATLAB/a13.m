%% A.13

x_sent = [XIn; XQn];
x_detected = [detected_XI; detected_XQ];

errors_num = 0;
for i = 1:N
    for j = 1:2
        difference(j,i) = x_sent(j,i) - x_detected(j,i);
        if (difference(j,i) ~= 0)
            errors_num = errors_num + 1;
        end
    end
end
errors_num