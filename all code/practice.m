% OFDM single transmitter/receiver transmission
% Name: Benjamin Pham
% Student ID: 25957066
% 
% Task: to simulate a transmission of data using OFDM technique
% blocks to so simulate, modulation block, IFFT, P/S,
% FFT, demodulate, recover signal. then add (2) Noise,(3) Multipath Channel

clc;clear all;close all
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
snr = [0:1:40];

% attenuation
% attenuation = 1;
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
t_data = round(rand(bits,1));                          % generate data
t_data = dec2bin(t_data);                              % changes vector to char type

%% symbol mapping
% 4QAM, 16QAM

reshape_data = reshape(t_data,length(t_data)/symbols, symbols );                   % reshape into pairings of 2 bits per symbol

dec_data = bin2dec(reshape_data);                                                  % must be a char vector

if mod_value == 1
    mod_data = qammod(dec_data,4,'unitaveragepower',true);
else 
    mod_data = qammod(dec_data,16,'unitaveragepower',true);
end 
%figure()    
%scatterplot(mod_data)

%% Modulation

% n = 0:pi/S:2*pi-pi/S;
% 
% in_phase = cos(n+pi/4);
% quadrature = sin(n+pi/4);
% symbol_book = (in_phase + quadrature*1i)';
% 
% X = symbol_book(dec_data+1);
% scatterplot(X)
%% IFFT
% Moves the data from the freq domain to the time domain

X = mod_data;   
X_blocks = reshape(X,n_fft,length(X)/n_fft);                        % reshape into 64 block subcarriers

x = ifft(X_blocks);                                                 % Inverse fast fourier transform

%% add CP
x_cp = [x(:,(end - n_cp + 1):end),x];                                     % add CP to start of data  

% x_cp2 = x_cp(:,(n_cp + 1):end);     
%% Parallel data to Serial stream

x_s = x_cp(:);

%%                          CHANNEL
%% Multipath Channel
% delays and attentuation that affects 
%channel = attenuation;      

%H_x = conv(x_s,channel);

H_x = x_s;
%% AWGN Noise
for i = 1:length(snr)
    H_noise = awgn(H_x,snr(i));

%%                          RECEIVER
% Serial to Parallel
% shapes the data back to 
y_p = reshape(H_noise,n_fft, length(H_noise + n_cp - 1)/n_fft);          

% remove cp
Y_rcp = y_p(:,(n_cp + 1):end);

%% FFT
% converts signal from time domain to frequency domain
Y_blocks = fft(Y_rcp);


%% Channel Estimation
% because we send pilot symbols, we can estimate the channel. symbols that
% are known beforehand. So its the received signal divided by the pilot
% symbol.

%H_hat = Y(1:n_p)./X(1:n_p);                 % estimating channel using pilot symbols

%Y_hat = H_hat .* X(1:n_p);                  


% mean squared error
% for m = 1:length(H_hat)
%     se(m) = (abs(Y(m) - Y_hat(m)))^2;
%     mse(m) = sum(se)/length(se);
% end
% mmse(mod_value,i) = min(mse);

%% equalisation???
% H_hate = abs(mean(H_hat));
% 
% Y_blocks2 = Y_blocks/H_hate;                % cancelling out the channel effects?

%% Demodulate
if mod_value == 1
    y = qamdemod(Y_blocks,4,'unitaveragepower',true);
else 
    y = qamdemod(Y_blocks,16,'unitaveragepower',true);
end

received_sym = dec2bin(y);
received_sig = reshape(received_sym, length(t_data), 1);


%% calculate BER
errors(mod_value,i) = 0;

for k = [1:length(received_sig)]
    if received_sig(k) ~= t_data(k)
        errors(mod_value,i) = errors(i) + 1;
    end
end

ber(mod_value,i) = errors(mod_value,i)/length(t_data);


end

% figure(mod_value)
semilogy(snr,ber(mod_value,:),'x-');
legend('4QAM','16QAM')
title('BER vs SNR'); 
xlabel('SNR'); ylabel('BER'); grid on; hold on

% subplot(2,2,2*mod_value)
% semilogy(eb_no,ber(mod_value,:));
% title('BER vs Eb/No'); legend('4QAM','16QAM')
% xlabel('Eb/No'); ylabel('BER'); grid on; hold on
end