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
n = 16;
bits = 2^n;
p = log2(bits);

% Modulation type (4QAM or 16QAM)
mod_type = '4QAM';

% bits per symbol
if mod_type == '4QAM'
    S = 2;   
else
    S = 4;
end

% fft/ifft size
n_fft = 64;

% cyclic prefix size
n_cp = 16;

% snr
snr = [0:1:30];

%%                          TRANSMITTER

%% Generate data to be sent

% 64kb of data 

t_data = round(rand(bits,1));                                       % generate data

t_data = dec2bin(t_data);

%% symbol mapping
% BPSK, 4QAM, 16QAM, 64 QAM

cons_data = reshape(t_data,length(t_data)/S, S );                   % reshape into pairings of 2 bits per symbol

dec_data = bin2dec(cons_data);                                      % must be a char vector

if mod_type == '4QAM'
    mod_data = qammod(dec_data,4,'unitaveragepower',true);
else 
    mod_data = qammod(dec_data,16,'unitaveragepower',true);
end 
    
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
% Moves the data from the time domain to the frequency domain

X = mod_data;   
X_blocks = reshape(X,n_fft,length(X)/n_fft);                        % reshape into 64 block subcarriers

x = ifft(X_blocks);                                                 % Inverse fast fourier transform

%% add CP
x_cp = [x(end-n_cp+1:end,:);x];                                     % add CP to end of data                                 
%% Parallel data to Serial stream

x_s = x_cp(:);

%%                          CHANNEL
%% Multipath Channel



%% AWGN Noise
for i = 1:length(snr)
    x_noise = awgn(x_s,snr(i),'measured');

%%                          RECEIVER
%% Serial to Parallel

x_p = reshape(x_noise,n_fft+n_cp, length(x_noise)/(n_fft+n_cp));

% remove cp
x_p_cp = x_p(n_cp + 1:end,:);

%% FFT

X_hat_blocks = fft(x_p_cp);

%% Demodulate
if mod_type == '4QAM'
    x_hat = qamdemod(X_hat_blocks,4,'unitaveragepower',true);
else 
    x_hat = qamdemod(X_hat_blocks,16,'unitaveragepower',true);
end

rec_sym = dec2bin(x_hat);

rec_sig = reshape(rec_sym, length(t_data), 1);

count(i) = 0;

for k = [1:length(rec_sig)]
    if rec_sig(k) ~= t_data(k)
        count(i) = count(i) + 1;
    end
end

ber(i) = count(i)/length(t_data);
end

semilogy(snr,ber);
title('BER vs SNR using 4QAM');
xlabel('SNR'); ylabel('BER'); grid on;