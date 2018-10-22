clc; clear all;close all; warning off
% simulate 2x1 OFDMA System with timing offsets 

data = 2^16;                                    % data points

n_fft = 128;                                     % fft size
    
n_cp = 16;                                      % cyclic prefix

snr = [0:2:30];

errors1 = zeros(size(snr));
ber1 = zeros(size(snr));
errors2 = zeros(size(snr));
ber2 = zeros(size(snr));

OFDM = n_fft+2*n_cp;
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
Xk1 = reshape(Tx1,n_fft/2,length(Tx1)/(n_fft/2)); % 64 sub carriers, per user
Xk2 = reshape(Tx2,n_fft/2,length(Tx2)/(n_fft/2));



%% ifft
% ofdm symbol

% ifft
    Xn1 = ifft(Xk1);
    Xn2 = ifft(Xk2);
% CP
    Xn1_cp = [Xn1((end - n_cp + 1):end,:);Xn1];  
    Xn2_cp = [Xn2((end - n_cp + 1):end,:);Xn2];  
  
    
Xn1_pad = [ Xn1_cp ; zeros(size(Xn1_cp)) ];
Xn2_pad = [ zeros( size(Xn2_cp)) ; Xn2_cp ];

% Parallel to Serial
    xn1 = Xn1_pad(:);
    xn2 = Xn2_pad(:);

%% channel


for k = 1:length(snr)
% channel gain
    yn11 = xn1*h(1);            % y time, receiver 
    yn21 = xn2*h(2);
db = snr(k);
    %noise
    yn11 = awgn(yn11,db,'measured');
    yn21 = awgn(yn21,db,'measured');
%    % delay
    delay11 = 23;
    delay21 = 0;
    
    
    % tail of previous symbol added to front of next symbol - ISI
    for i = 0:(length(yn11)/(OFDM))-2
        if delay11 ~= 0
        yn11(1+OFDM+OFDM*i:1+OFDM+delay11+OFDM*i) = yn11(1+OFDM+OFDM*i:1+OFDM+delay11+OFDM*i) ...
            + yn11(OFDM/2-delay11+OFDM*i:OFDM/2+OFDM*i);
        end
        if delay21 ~= 0
        yn21(1+OFDM+OFDM*i:1+OFDM+delay21+OFDM*i) = yn21(1+OFDM+OFDM*i:1+OFDM+delay21+OFDM*i) ...
            + yn21(OFDM-delay21+OFDM*i:OFDM+OFDM*i);
        end
    
    end
    if delay11 ~= 0
    yn11(1:OFDM) = [zeros(delay11+1,1);yn11(1:OFDM-delay11-1)];      
    end
    if delay21 ~= 0
    yn21(OFDM+1:2*OFDM) = [zeros(delay21+1,1);yn21(OFDM+1:2*OFDM-delay21-1)];  
    end
%         superposition
    yn = yn11 + yn21; 
    
%% receiver
    %serial to parallel 128+CP streams
    yn_sub = reshape(yn,n_fft+2*n_cp,length(Xn1));
% %   users
    Yn_sub1 = yn_sub(1:80,:);                               % higher order decision block
    Yn_sub2 = yn_sub(81:160,:);
    
%     % delay
%     delay11 = 20;
%     delay21 = 0;
% 
%     for i = 0:(length(yn11)/(OFDM))-2
%         if delay11 ~= 0
%         Yn11(1+OFDM+OFDM*i:1+OFDM+delay11+OFDM*i) = Yn11(1+OFDM+OFDM*i:1+OFDM+delay11+OFDM*i) ...
%             + yn11(OFDM-delay11+OFDM*i:OFDM+OFDM*i);
%         end
%         if delay21 ~= 0
%         yn21(1+OFDM+OFDM*i:1+OFDM+delay21+OFDM*i) = yn21(1+OFDM+OFDM*i:1+OFDM+delay21+OFDM*i) ...
%             + yn21(OFDM-delay21+OFDM*i:OFDM+OFDM*i);
%         end
% 
%     end
%     if delay11 ~= 0
%     yn11(1:OFDM) = [zeros(delay11+1,1);yn11(1:OFDM-delay11-1)];      
%     end
%     if delay21 ~= 0
%     yn21(1:OFDM) = [zeros(delay21+1,1);yn21(1:OFDM-delay21-1)];  
%     end    
    
        % remove CP
    yn_sub_rcp1 = Yn_sub1((n_cp+1):end,:);
    yn_sub_rcp2 = Yn_sub2((n_cp+1):end,:);
    
    %% DFT convert back to freq domain
    Yk_block1 = fft(yn_sub_rcp1);
    Yk_block2 = fft(yn_sub_rcp2);

% P/S
    Yk1_s = Yk_block1(:);              % transform 64 subchannels back into 1 stream
    Yk2_s = Yk_block2(:);       
    
%     delay = 0;
%     Yk1_s = [zeros(delay,1); Yk1_s];
    
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


errors1(mod_value,k) = 0;
for i = 1:length(binary_data1)-n_fft*symbols
    if output1(i+n_fft*symbols) ~= binary_data1(i+n_fft*symbols)
        errors1(mod_value,k) = errors1(mod_value,k) + 1;
    end
end
ber1(mod_value,k) = errors1(mod_value,k)/length(binary_data1);


errors2(mod_value,k) = 0;
for i = 1:length(binary_data2)-n_fft*symbols
    if output2(i+n_fft*symbols) ~= binary_data2(i+n_fft*symbols)
        errors2(mod_value,k) = errors2(mod_value,k) + 1;
    end
end
ber2(mod_value,k) = errors2(mod_value,k)/length(binary_data2);


end
semilogy(snr-10*log10(symbols),ber1(mod_value,:),'-',snr-10*log10(symbols),ber2(mod_value,:),'-')
%semilogy(snr,ber2(mod_value,:),'o-');
title('Two users One receiver with Time Offset greater than CP by 50%'); legend('4QAM user1','4QAM user2','16QAM user1','16QAM user2','64QAM user1','64QAM user2');
xlabel('E_b/N_o (dB)'); ylabel('BER'); grid on;hold on
end
error = zeros(1,length(binary_data1));
for i = 1:length(binary_data1)
    if output1(i) ~= binary_data1(i)
        error(i) = 1;
    end
end

% scatterplot(Yk_block1(:,1))
% 
% scatterplot(Yk_block1(:,2))