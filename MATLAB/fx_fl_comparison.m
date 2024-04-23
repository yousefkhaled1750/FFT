
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
W_D_NoOfBits = 10;  
W_C_NoOfBits = 8;  
for W_C_NoOfBits = 6:12
    for W_D_NoOfBits = 8:16
        W_fx = floor(W*2^W_C_NoOfBits)/2^W_C_NoOfBits;
        previous_stage_fl = x_n_reorder;
        previous_stage_fx = x_n_reorder *2^(W_D_NoOfBits-numOfBits_in); %shift left to maintain the scaling of W_D_NoOfBits
        current_stage_fl  = zeros(1,N);
        current_stage_fx  = zeros(1,N);
        for s = 1:log2(N)
        %performing the stage operations
        for i = 0:N/2^s - 1
            for j = 1 : 2^(s-1)
                % the multiplication part
                temp_fl = previous_stage_fl(2^s*i + 2^(s-1) + j)*W((j-1)*N/2^s + 1);            
                previous_stage_fl(2^s*i + 2^(s-1) + j) = temp_fl;
    
                temp_fx = previous_stage_fx(2^s*i + 2^(s-1) + j)*W_fx((j-1)*N/2^s + 1);
                % quantizing the output into the potential wordlength 
                temp_fx = floor(temp_fx*2^W_C_NoOfBits)/2^W_C_NoOfBits;
                previous_stage_fx(2^s*i + 2^(s-1) + j) = temp_fx;
    
                current_stage_fl(2^s*i + j)           =  previous_stage_fl(2^s*i + j) + previous_stage_fl(2^s*i+j + 2^(s-1));
                current_stage_fl(2^s*i + j + 2^(s-1)) =  previous_stage_fl(2^s*i + j) - previous_stage_fl(2^s*i+j + 2^(s-1));
                
                op1 = previous_stage_fx(2^s*i + j);
                op2 = previous_stage_fx(2^s*i+j + 2^(s-1));
                temp1           =  op1 + op2;
                temp2           =  op1 - op2;
                % eliminate the additional bit from lsb
                temp1 = floor(temp1/2)*2;
                temp2 = floor(temp2/2)*2;
                % shift right by the W_D fractional bits
                temp1 = temp1/2^(W_D_NoOfBits-numOfBits_in);
                temp2 = temp2/2^(W_D_NoOfBits-numOfBits_in);
                current_stage_fx(2^s*i + j) = temp1;
                current_stage_fx(2^s*i + j + 2^(s-1)) = temp2;
    
            end
        end
        previous_stage_fl = current_stage_fl;
        previous_stage_fx = current_stage_fx*2^(W_D_NoOfBits-numOfBits_in);
        end
        % here we scale down the registers to the actual values scaling
        X_f_Qunatized = previous_stage_fx/2^(W_D_NoOfBits-numOfBits_in);
        X_f_fl= previous_stage_fl;
        
        X_f_fl_power = abs(X_f_fl).^2;
        X_f_fl_power_Avg = sum(X_f_fl_power)/N;
        fx_error_power = abs(X_f_fl - X_f_Qunatized).^2;
        fx_error_power_avg = sum(fx_error_power)/N;
        SQNR(W_C_NoOfBits,W_D_NoOfBits) = pow2db(X_f_fl_power_Avg/fx_error_power_avg);
    end
end
% here we scale down the registers to the actual values scaling
% X_f_Qunatized = previous_stage_fx/2^(W_D_NoOfBits-numOfBits_in);
% X_f_fl= previous_stage_fl;
% 
% X_f_fl_power = abs(X_f_fl).^2;
% X_f_fl_power_Avg = sum(X_f_fl_power)/N;
% fx_error_power = abs(X_f_fl - X_f_Qunatized).^2;
% fx_error_power_avg = sum(fx_error_power)/N;
% SQNR = pow2db(X_f_fl_power_Avg/fx_error_power_avg)

x = linespace(8,1,16);
y = linespace(6,1,12);
mesh(SQNR,)