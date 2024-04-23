% Dummy implementation of an N-point DIT FFT algorithm using butterfly flow

clear ;
close all;
clc;
N = 1024;
n = (0:1:N-1);
x_n = [ones([1,N/4]), -1*ones([1,N/4]), ones([1,N/4]), -1*ones([1,N/4])];
x_n_reorder = x_n;
%construct the twiddle factor 
W = exp(-1i*2*pi*n(1:N/2)/N);   % we only need half of the twiddles


%% applying the bit reverse algorithm
i = 1;
for r = 1:N
    if r > i
        % swap x_n and x_i
        temp = x_n_reorder(r);
        x_n_reorder(r) = x_n_reorder(i);
        x_n_reorder(i) = temp;
    end
    m = N/2;
    while m >= 2 && i > m
        i = i - m;
        m = m/2;
    end
    i = i+m;
end

%% performing the subsequent stages
previous_stage = x_n_reorder;
current_stage  = zeros(1,N);
for s = 1:log2(N)
    %performing the stage operations
    for i = 0:N/2^s - 1
        for j = 1 : 2^(s-1)
            previous_stage(2^s*i + 2^(s-1) + j) = previous_stage(2^s*i + 2^(s-1) + j)*W((j-1)*N/2^s + 1);
            current_stage(2^s*i + j)           =  previous_stage(2^s*i + j) + previous_stage(2^s*i+j + 2^(s-1));
            current_stage(2^s*i + j + 2^(s-1)) =  previous_stage(2^s*i + j) - previous_stage(2^s*i+j + 2^(s-1));
        end
    end
    previous_stage = current_stage;
end
X_f = current_stage;
X_f_mat = fft(x_n,N);

figure;
plot(real(X_f))
hold on
plot(imag(X_f))
legend('real','imag')
title("Danielson-Lanczos algorithm")

figure
plot(real(X_f_mat))
hold on
plot(imag(X_f_mat))
legend('real','imag')
title("MATLAB implementation")

