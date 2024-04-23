% the specs: 64-pt fft radix 2
%            SQNR = 40 dB
%            input s(8,0)   signed pure integer
% the problem is to choose the size of the twiddle factors and the growth
% bits in the subsequent stages to get the desired SQNR
% the methodology is to try out different word lengths W_D in FFT with constant
% twiddle factors size W_C and then try different twiddle factors sizes 
% we can use nested loops (the outer is W_C and the inner is W_D)

clear ;
close all;
clc;
numOfBits_in = 8;
N = 64;
n = (0:1:N-1);
x_n = rand([1 N]);

% quantization of input
x_n = floor((2*x_n - 1)*((2^numOfBits_in) - 1));

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
W_D_NoOfBits = 12;  % it would be 8 bits as the input and 4 bits fractional
                    % and when we iterate over W_D we will modify only the fractional bits
% we need first to quantize the twiddle factors
% starting with 8 bits as the input wordlength 
W_C_NoOfBits = 7;   % it's actually 8 but we assigend 7 fractional bits and 1 integral bits for the "1.000"
W_fx = floor(W*2^W_C_NoOfBits)/2^W_C_NoOfBits;
previous_stage = x_n_reorder *2^(W_D_NoOfBits-numOfBits_in); %shift left to maintain the scaling of W_D_NoOfBits
current_stage  = zeros(1,N);
for s = 1:log2(N)
    %performing the stage operations
    for i = 0:N/2^s - 1
        for j = 1 : 2^(s-1)
            % the multiplication part
            % previous_stage(2^s*i + 2^(s-1) + j) = previous_stage(2^s*i + 2^(s-1) + j)*W((j-1)*N/2^s + 1);            
            temp = previous_stage(2^s*i + 2^(s-1) + j)*W_fx((j-1)*N/2^s + 1);
            % the previous operation caused the integral bits to be 8+1 and
            % the fractionals bits to be W_C + (W_D - 8)
            % so we need to quantize it back to W_D
            % by removing the W_C expansion

            % quantizing the output into the potential wordlength 
            temp = floor(temp*2^W_C_NoOfBits)/2^W_C_NoOfBits;
            previous_stage(2^s*i + 2^(s-1) + j) = temp;
            % to add 2 W_D (8+x) bits registers, we need first to shift
            % left both by x, add and then shift right by 1 and (floor) where we
            % eliminate the additional bit of the operation from the lsb
            % then shift by x
            op1 = previous_stage(2^s*i + j);
            op2 = previous_stage(2^s*i+j + 2^(s-1))*2^(W_D_NoOfBits-numOfBits_in);

            % current_stage(2^s*i + j)           =  previous_stage(2^s*i + j) + previous_stage(2^s*i+j + 2^(s-1));
            % current_stage(2^s*i + j + 2^(s-1)) =  previous_stage(2^s*i + j) - previous_stage(2^s*i+j + 2^(s-1));
            temp1           =  op1 + op2;
            temp2           =  op1 - op2;
            % eliminate the additional bit from lsb
            temp1 = floor(temp1/2);
            temp2 = floor(temp2/2);
            % shift right by the W_D fractional bits
            temp1 = temp1/2^(W_D_NoOfBits-numOfBits_in);
            temp2 = temp2/2^(W_D_NoOfBits-numOfBits_in);
            current_stage(2^s*i + j) = temp1;
            current_stage(2^s*i + j + 2^(s-1)) = temp2;

        end
    end
    previous_stage = current_stage*2^(W_D_NoOfBits-numOfBits_in);
end
% here we scale down the registers to the actual values scaling
X_f_Qunatized = current_stage/2^(W_D_NoOfBits-numOfBits_in);
X_f_mat = fft(x_n,N);

X_f_mat_power = abs(X_f_mat).^2;
X_f_mat_power_Avg = sum(X_f_mat_power)/N;
fx_error_power = abs(X_f_mat - X_f_Qunatized).^2;
fx_error_power_avg = sum(fx_error_power)/N;
SQNR = pow2db(X_f_mat_power_Avg/fx_error_power_avg)

