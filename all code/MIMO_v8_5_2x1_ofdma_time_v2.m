clc; clear all;close all; warning off
% simulate 2x1 OFDMA System with timing offsets 

data = 2^16;                                    % data points

n_fft = 128;                                     % fft size
    
n_cp = 16;                                      % cyclic prefix

snr = [0:1:35];

%delay of unit length
delay = 2;

errors1 = zeros(size(snr));
ber1 = zeros(size(snr));
errors2 = zeros(size(snr));
ber2 = zeros(size(snr));
%% generate data

binary_data1 = round(rand(data,1));              % generate data
binary_data2 = round(rand(data,1));
bin1= binary_data1;
bin2= binary_data2;
% channel
h = [randn() + randn()*1i];
h = [h h];
%h = [randn()+randn()*1i , randn()+randn()*1i]; 
%% modulation

for mod_value = 1:3
    % bits per symbol
    if mod_value == 1
        symbols = 2;   
    elseif mod_value == 2
        symbols = 4;
    else 
        symbols = 6;
        while floor(length(binary_data1)/symbols) ~= length(binary_data1)/symbols
        binary_data1 = [binary_data1; zeros(1,1)];                                    %padding for symbol mapping
        end
        while floor(length(binary_data1)/symbols/n_fft) ~= length(binary_data1)/symbols/n_fft
        binary_data1 = [binary_data1; zeros(6,1)];                                    %padding for subcarriers mapping
        end
    
        while floor(length(binary_data2)/symbols) ~= length(binary_data2)/symbols
        binary_data2 = [binary_data2; zeros(1,1)];                                    %padding for symbol mapping
        end
        while floor(length(binary_data2)/symbols/n_fft) ~= length(binary_data2)/symbols/n_fft
        binary_data2 = [binary_data2; zeros(6,1)];                                    %padding for subcarriers mapping
        end
    
    end
    
mod_method = 2^symbols;       
                                                
mod_data1 = qammod(binary_data1,mod_method,'unitaveragepower',true,'inputtype','bit');
mod_data2 = qammod(binary_data2,mod_method,'unitaveragepower',true,'inputtype','bit');

%% splitting data up between the transmitters
Tx1 = mod_data1;
Tx2 = mod_data2;

% allocating data onto each subcarrier

%% subcarriers
Xk1 = reshape(Tx1,n_fft/2,length(Tx1)/(n_fft/2)); % 64 sub carriers, 
Xk2 = reshape(Tx2,n_fft/2,length(Tx2)/(n_fft/2));

Xk1_pad = [ Xk1 ; zeros(size(Xk1)) ];
Xk2_pad = [ zeros( size(Xk2)) ; Xk2 ];

%% ifft
% ofdm symbol
for s = 1:length(Xk1_pad)
% ifft
    Xn1(:,s) = ifft(Xk1_pad(:,s));
    Xn2(:,s) = ifft(Xk2_pad(:,s));
% CP
    Xn1_cp(:,s) = [Xn1((end - n_cp + 1):end,s);Xn1(:,s)];  
    Xn2_cp(:,s) = [Xn2((end - n_cp + 1):end,s);Xn2(:,s)];  
  
end
% Parallel to Serial
    xn1 = Xn1_cp(:);
    xn2 = Xn2_cp(:);

%% channel


for k = 1:length(snr)
% channel gain
    yn11 = xn1*h(1);            % y time, receiver 
    yn21 = xn2*h(2);

    %noise
    yn11 = awgn(yn11,snr(k),'measured');
    yn21 = awgn(yn21,snr(k),'measured');
% 
%     % time delay
%     yn11 = [zeros(delay,1) ; yn11 ];
%     yn21 = [zeros(delay,1) ; yn21 ];
    
    % superposition
    yn = yn11 + yn21; 
    
%     yn2 = yn(1:(end-delay));
    %serial to parallel 128 streams
    yn_sub = reshape(yn,n_fft+n_cp,length(Xn1));
    % remove CP
    yn_sub_rcp = yn_sub((n_cp+1):end,:);
    
    %% DFT convert back to freq domain
    Yk_block = fft(yn_sub_rcp);

    % decode sub bands to seperate users
    Yk1 = Yk_block(1:64,:);                               % higher order decision block
    Yk2 = Yk_block(65:128,:);                             % shifting the data, causes the information to go 
                                                          % into the wrong
                                                          % subbands
    Yk1_s = Yk1(:);              % transform 64 subchannels back into 1 stream
    Yk2_s = Yk2(:);            
    Yk1_s = [zeros(delay,1); Yk1_s];
    Yk2_s = [zeros(delay,1); Yk2_s];
%% Channel estimation
% assume perfect channel estimation, so we know what h1 and h2 are

Xhat1 = Yk1_s./h(1);
Xhat2 = Yk2_s./h(2);

%% demodulate
X_demod1 = qamdemod(Xhat1,mod_method,'unitaveragepower',true,'outputtype','bit');
X_demod2 = qamdemod(Xhat2,mod_method,'unitaveragepower',true,'outputtype','bit');

output1 = X_demod1(:).';
output2 = X_demod2(:).';
%%
%  binary_data1 = [bin1; zeros(delay,1)];
%  binary_data2 = [bin2; zeros(delay,1)];

output1 = [output1, zeros(1,symbols*delay)];
output2 = [output2, zeros(1,symbols*delay)];

errors1(mod_value,k) = 0;
for i = 1:length(binary_data1)
    if output1(i+symbols*delay) ~= binary_data1(i)
        errors1(mod_value,k) = errors1(mod_value,k) + 1;
    end
end
ber1(mod_value,k) = errors1(mod_value,k)/length(binary_data1);


errors2(mod_value,k) = 0;
for i = 1:length(binary_data2)
    if output2(i+symbols*delay) ~= binary_data2(i)
        errors2(mod_value,k) = errors2(mod_value,k) + 1;
    end
end
ber2(mod_value,k) = errors2(mod_value,k)/length(binary_data2);


end
semilogy(snr,ber1(mod_value,:),'x-',snr,ber2(mod_value,:),'o-')
%semilogy(snr,ber2(mod_value,:),'o-');
title('OFDMA for 2x1 system'); legend('4QAM user1','4QAM user2','16QAM user1','16QAM user2','64QAM user1','64QAM user2');
xlabel('SNR'); ylabel('BER'); grid on;hold on
end
error = zeros(1,length(binary_data1));
for i = 1:length(binary_data1)
    if output1(i) ~= binary_data1(i)
        error(i) = 1;
    end
end