% OFDM 2 transmitter/receiver transmission
% Name: Benjamin Pham
% Student ID: 25957066
% 
% Task: to simulate a transmission of data using OFDM technique
% blocks to so simulate, modulation block, IFFT, P/S,
% FFT, demodulate, recover signal. then add (2) Noise,(3) Multipath Channel

clc;clear all; close all; warning off
%% Set Parameters

% amount of data to be transmitted
% 64kb is 2^16
% ^20 is a megabit
n = 16;
bits = 2^n;
p = log2(bits);


% fft/ifft size
n_fft = 64;

% cyclic prefix size
n_cp = 0.25*n_fft;

% snr
snr = [0:1:30];

% attenuationx

% eb_no = snr.*(B/Rb);
t_data = round(rand(bits,1));                                              % generate data
     
OFDM = n_fft+n_cp;
channel = randn() +j*randn()

DS = 63;
PS = 1;
% NS = n_fft - DS - PS;

for mod_value = [1:3]
% bits per symbol
    if mod_value == 1
        symbols = 2; 
        pilot = [0 1]';
                while floor(length(t_data)/symbols) ~= length(t_data)/symbols
                t_data = [t_data; zeros(1,1)];                                    %padding for symbol mapping
                end
                while floor(length(t_data)/symbols/DS) ~= length(t_data)/symbols/DS
                t_data = [t_data; zeros(1,1)];                                    %padding for subcarriers mapping
                end                
    elseif mod_value == 2
        symbols = 4;
        pilot = [0 0 0 1]';
                while floor(length(t_data)/symbols) ~= length(t_data)/symbols
                t_data = [t_data; zeros(1,1)];                                    %padding for symbol mapping
                end
                while floor(length(t_data)/symbols/DS) ~= length(t_data)/symbols/DS
                t_data = [t_data; zeros(1,1)];                                    %padding for subcarriers mapping
                end                  
    elseif mod_value == 3
        symbols = 6;
        pilot = [0 0 0 0 0 1]';
                while floor(length(t_data)/symbols) ~= length(t_data)/symbols
                t_data = [t_data; zeros(1,1)];                                    %padding for symbol mapping
                end
                while floor(length(t_data)/symbols/DS) ~= length(t_data)/symbols/DS
                t_data = [t_data; zeros(1,1)];                                    %padding for subcarriers mapping
                end
    
    end
    
 mod_method = 2^symbols;        
%%                          TRANSMITTER
%% Generate data to be sent

% 64kb of data 
                                             % changes vector to char type

%% symbol mapping
% 4QAM, 16QAM


mod_data = qammod(t_data,mod_method,'unitaveragepower',true,'inputtype','bit');
pilot = qammod(pilot,mod_method,'unitaveragepower',true,'inputtype','bit');

%figure()    
%scatterplot(mod_data)

%% IFFT
% Moves the data from the time domain to the frequency domain

X = mod_data;   

%         while floor(length(X)/DS) ~= length(X)/DS
%         X = [X; zeros(1,1)];                                    %padding for symbol mapping
%         end

% data symbols in OFDM frame
X_data = reshape(X,DS,length(X)/DS);            
% pilot insertion
X_blocks = [ones(1,length(X_data))*pilot;X_data];

% IDFT operation
x = ifft(X_blocks);                                                 % Inverse fast fourier transform

%% add CP
x_cp = [x((end - n_cp + 1):end,:);x];                                     % add CP to end of data                                 
%% Parallel data to Serial stream

x_s = x_cp(:);

%%                          CHANNEL
%% Multipath Channel
% delays and attentuation that affects 

H_x = x_s*channel;

%H_x = x_s;

%% AWGN Noise
for i = 1:length(snr)
    H_noise = awgn(H_x,snr(i),'measured');

    
    delay = 17;
    % tail of previous symbol added to front of next symbol - ISI
    for j = 0:(length( H_noise)/(OFDM))-2
        if delay ~= 0
         H_noise(1+OFDM+OFDM*j:1+OFDM+delay+OFDM*j) =  H_noise(1+OFDM+OFDM*j:1+OFDM+delay+OFDM*j) ...
            +  H_noise(OFDM-delay+OFDM*j:OFDM+OFDM*j);
        end
    end
    
    if delay ~= 0
     H_noise(1:OFDM) = [zeros(delay+1,1); H_noise(1:OFDM-delay-1)];      
    end

    
%%                          RECEIVER
%% Serial to Parallel
y_p = reshape(H_noise, OFDM, length(H_noise)/OFDM);

% remove cp
y_p_cp = y_p((n_cp + 1):end,:);

%% FFT
% converts signal from time domain to frequency domain
Y_blocks = fft(y_p_cp);


Y_pilots = Y_blocks(1,:);
% extract data symbols
Y_data = Y_blocks(2:64,:);

Y = Y_data(:);
%% Channel Estimation
% because we send pilot symbols, we can estimate the channel. symbols that
% are known beforehand. So its the received signal divided by the pilot
% symbol.

H_hat = Y_pilots./pilot;                 % estimating channel using pilot symbols



%% equalisation
H_hate = mean(H_hat);                       % found a mistake where taking the abs of it caused a very weird results in the equalisation

Y_blocks2 = Y./H_hate;                % cancelling out the channel effects
                                             

%% Demodulate equalistation part

y2 = qamdemod(Y_blocks2,mod_method,'unitaveragepower',true,'outputtype','bit');

received_sig2 = y2(:);
errors2(mod_value,i) = 0;
% received_sig2b = [zeros(delay,1); received_sig2];


for m = [1:length(received_sig2)-n_fft*symbols]
    if received_sig2(m+n_fft*symbols) ~= t_data(m+n_fft*symbols)
        errors2(mod_value,i) = errors2(mod_value,i) + 1;
    end
end

ber2(mod_value,i) = errors2(mod_value,i)/length(t_data);
end
%figure()
%semilogy(snr,ber(mod_value,:),'x-',snr,ber2(mod_value,:),'o-');
semilogy(snr-10*log10(symbols),ber2(mod_value,:),'-');
title('BER vs SNR'); legend('4QAM','16QAM','64QAM');
xlabel('E_b/N_o (dB)'); ylabel('BER'); grid on; hold on


end
