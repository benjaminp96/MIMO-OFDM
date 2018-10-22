% OFDM 2 transmitter/receiver transmission
% Name: Benjamin Pham
% Student ID: 25957066
% 
% Task: to simulate a transmission of data using OFDM technique
% blocks to so simulate, modulation block, IFFT, P/S,
% FFT, demodulate, recover signal. then add (2) Noise,(3) Multipath Channel

clc;clear all;close all; warning off
%% Set Parameters

% amount of data to be transmitted
% 64kb is 2^16
% ^20 is a megabit
n = 16;
bits = 2^n;
p = log2(bits);



% Modulation type (4QAM or 16QAM)           4QAM = 1;  16QAM = 2;
mod_type = '4QAM';
mod_types = {'4QAM','16QAM'};
mod_value = find (ismember (mod_types, mod_type));

% Eb/No             assume 4G network, SNR = Eb/No * Rb/B           
                    % rb is bit rate = 100Mbps , B = 20MHz
% Rb = 100*10^6;
% B = 20*10^6;

% fft/ifft size
n_fft = 64;
%pilot length
n_p = 4;

% cyclic prefix size
n_cp = 16;

% snr
snr = [0:1:40];

% attenuation
attenuation = rand() + j*rand();
% eb_no = snr.*(B/Rb);

for mod_value = [1:2]
% bits per symbol
if mod_value == 1
    symbols = 2;   
else
    symbols = 4;
end
%%                          TRANSMITTER
%% Generate data to be sent

% 64kb of data 
t_data = round(rand(bits,1));                                              % generate data
t_data = dec2bin(t_data);                                                  % changes vector to char type

%% symbol mapping
% 4QAM, 16QAM

reshape_data = reshape(t_data,length(t_data)/symbols, symbols );                       % reshape into pairings of 2 bits per symbol

decimal_data = bin2dec(reshape_data);                                      % must be a char vector

if mod_value == 1
    mod_data = qammod(decimal_data,4,'unitaveragepower',true);
else 
    mod_data = qammod(decimal_data,16,'unitaveragepower',true);
end 
%figure()    
%scatterplot(mod_data)


%%
X = mod_data;   
X_block = reshape(X,n_fft,length(X)/n_fft);                        % reshape into 64 block subcarriers

%% insert pilot symbols
pilots = ones(1,1);

if mod_value == 1
    mod_pilots = qammod(pilots,4,'unitaveragepower',true);
else 
    mod_pilots = qammod(pilots,16,'unitaveragepower',true);
end 
[c,d] = size(X_block);

X_blocks = zeros(c,d+1);
for i2 = 1:n_fft
X_blocks(i2,:) = [mod_pilots(1,1), X_block(i2,:)];  
end 
% mod_pilots(1,2) X_block(i2,(length(X_block)/4+1):length(X_block)/2) ...
% mod_pilots(1,3) X_block(i2,(length(X_block)/2+1):3*(length(X_block)/4)) mod_pilots(1,4) X_block(i2,(3*(length(X_block)/4)+1):end)];

%% IFFT
% Moves the data from the time domain to the frequency domain
x = ifft(X_blocks);                                                 % Inverse fast fourier transform

%% add CP
x_cp = [x(:,(end - n_cp + 1) :end),x];                                     % add CP to end of data                                 
%% Parallel data to Serial stream

x_s = x_cp(:);

%%                          CHANNEL
%% Multipath Channel
% delays and attentuation that affects 
channel = attenuation;      

H_x = conv(x_s,channel);
%H_x = x_s;

%% AWGN Noise
for i = 1:length(snr)
    H_noise = awgn(H_x,snr(i),'measured');

%%                          RECEIVER
%% Serial to Parallel
y_p = reshape(H_noise, n_fft, length(H_noise + n_cp - 1)/n_fft);

% remove cp
x_p_cp = y_p(:,(n_cp + 1):end);

%% FFT
% converts signal from time domain to frequency domain
Y_blocks = fft(x_p_cp);
Y = reshape(Y_blocks, (length(X)+n_fft*length(pilots)),1);
%% Channel Estimation
% because we send pilot symbols, we can estimate the channel. symbols that
% are known beforehand. So its the received signal divided by the pilot
% symbol.

Y_hat = zeros(c,d+1);
for i2 = 1:n_fft
H_hat(i2) = Y_blocks(i2,1)./X_blocks(i2,1);                                     % estimating channel using pilot symbols

Y_hat(i2,:) = H_hat(i2) .* X_blocks(i2,:);                  
end

% % mean squared error
% for m = 1:length(H_hat)
%     se(m) = (abs(Y(m) - Y_hat(m)))^2;
%     mse(m) = sum(se)/length(se);
% end
% mmse(mod_value,i) = min(mse);

%% equalisation
H_hate = mean(H_hat);

Y_blocks2 = Y_blocks/H_hate;                % cancelling out the channel effects?

% remove pilot
Y_blocks = Y_blocks(:,2:end);
Y_blocks2 = Y_blocks2(:,2:end);

%% Demodulate
if mod_value == 1
    y = qamdemod(Y_blocks,4,'unitaveragepower',true);
else 
    y = qamdemod(Y_blocks,16,'unitaveragepower',true);
end

received_sym = dec2bin(y);

received_sig = reshape(received_sym, length(t_data), 1);

errors(mod_value,i) = 0;

for k = [1:length(received_sig)]
    if received_sig(k) ~= t_data(k)
        errors(mod_value,i) = errors(mod_value,i) + 1;
    end
end

ber(mod_value,i) = errors(mod_value,i)/length(t_data);

%% Demodulate equalistation part
if mod_value == 1
    y2 = qamdemod(Y_blocks2,4,'unitaveragepower',true);
else 
    y2 = qamdemod(Y_blocks2,16,'unitaveragepower',true);
end

received_sym2 = dec2bin(y2);
received_sig2 = reshape(received_sym2, length(t_data), 1);

errors2(mod_value,i) = 0;

for m = [1:length(received_sig2)]
    if received_sig2(m) ~= t_data(m)
        errors2(mod_value,i) = errors2(mod_value,i) + 1;
    end
end

ber2(mod_value,i) = errors2(mod_value,i)/length(t_data);
end
%figure()
semilogy(snr,ber(mod_value,:),'x-',snr,ber2(mod_value,:),'o-');
title('BER vs SNR'); legend('4QAM','4QAMeq','16QAM','16QAMeq');
xlabel('SNR'); ylabel('BER'); grid on; hold on


end
