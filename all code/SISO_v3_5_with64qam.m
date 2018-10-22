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

%pilot length
n_p = 0.1*bits;

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

% cyclic prefix size
n_cp = 16;

% snr
snr = [0:1:30];

% attenuationx
attenuation = rand()+j*rand();
attenuation = attenuation/abs(attenuation);
% eb_no = snr.*(B/Rb);

t_data = round(rand(bits,1));  

for mod_value = [1:3]
% bits per symbol
if mod_value == 1
    symbols = 2;   
elseif mod_value == 2
    symbols = 4;
else 
    symbols = 6;
    while floor(length(t_data)/symbols) ~= length(t_data)/symbols
    t_data = [t_data; zeros(1,1)];                                    %padding for symbol mapping
    end
    while floor(length(t_data)/symbols/n_fft) ~= length(t_data)/symbols/n_fft
    t_data = [t_data; zeros(6,1)];                                    %padding for subcarriers mapping
    end
    
end
%%                          TRANSMITTER

%% symbol mapping
% 4QAM, 16QAM

                                     % must be a char vector

mod_method = 2^symbols; 

mod_data = qammod(t_data,mod_method,'unitaveragepower',true,'inputtype','bit');
%figure()    
%scatterplot(mod_data)

%% IFFT
% Moves the data from the time domain to the frequency domain

X = mod_data;   
X_blocks = reshape(X,n_fft,length(X)/n_fft);                        % reshape into 64 block subcarriers

x = ifft(X_blocks);                                                 % Inverse fast fourier transform

%% add CP
x_cp = [x(:,(end - n_cp + 1) :end),x];                                     % add CP to end of data                                 
%% Parallel data to Serial stream

x_s = x_cp(:);

%%                          CHANNEL
%% Multipath Channel
% delays and attentuation that affects 
channel = 1; %attenuation;      

H_x = conv(x_s,channel,'same');
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
Y_blocks = fft(x_p_cp,n_fft);
Y = reshape(Y_blocks, length(X),1);
%% Channel Estimation
% because we send pilot symbols, we can estimate the channel. symbols that
% are known beforehand. So its the received signal divided by the pilot
% symbol.

H_hat = Y(1:n_p)./X(1:n_p);                 % estimating channel using pilot symbols

%Y_hat = H_hat .* X(1:n_p);                  


% % mean squared error
% for m = 1:length(H_hat)
%     se(m) = (abs(Y(m) - Y_hat(m)))^2;
%     mse(m) = sum(se)/length(se);
% end
% mmse(mod_value,i) = min(mse);

%% equalisation
H_hate = mean(H_hat);                       % found a mistake where taking the abs of it caused a very weird results in the equalisation

Y_blocks2 = Y_blocks./H_hate;                % cancelling out the channel effects
                                             
%% Demodulate

y = qamdemod(Y_blocks,mod_method,'unitaveragepower',true,'outputtype','bit');


received_sym = y;

received_sig = reshape(received_sym, length(t_data), 1);

errors(mod_value,i) = 0;

for k = [1:length(received_sig)]
    if received_sig(k) ~= t_data(k)
        errors(mod_value,i) = errors(mod_value,i) + 1;
    end
end

ber(mod_value,i) = errors(mod_value,i)/length(t_data);

%% Demodulate equalistation part
y2 = qamdemod(Y_blocks2,mod_method,'unitaveragepower',true,'outputtype','bit');


received_sym2 = y2;
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
%semilogy(snr,ber(mod_value,:),'x-',snr,ber2(mod_value,:),'o-');
semilogy(snr,ber2(mod_value,:),'o-');
title('Single Transmitter and Receiver System'); legend('4QAM','16QAM','64QAM');
xlabel('SNR'); ylabel('BER'); grid on; hold on


end
