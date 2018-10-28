%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Student:      Benjamin Pham
%   ID:           25957066
%   
%   Version:        Final
%   Changes:        
%   Description:    a 1x1 SISO OFDM System
%                                   
%   Simulates a 1x1 SISO OFDM system. This system includes the
%   Transmitter: 
%   QAM Modulation, S/P converter, Pilot Insertion IFFT, Add Cyclic Prefix, 
%   P/S converter
%   Channel: 
%   Frequency flat fading Channel modelled by Complex number, AWGN, 
%   Time Offset
%   Receiver: 
%   S/P converter, Removal of Cyclic Prefix, DFT, Channel estimation,  
%   Channel equalisation, P/S converter, QAM Demodulation
%
%   Simulation performs a BER test vs SNR when various time offsets are
%   inputted by sending binary data through the system and measuring the
%   amount of errors at the receiver. Time offsets can be changed beginning
%   on line 138 for different channels.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all; warning off

%% Input parameters
% data points
data = 2^16;                                   

% fft size
n_fft = 64;                                     

% cyclic prefix length
n_cp = 0.25*n_fft;                                  % CP length 25% of FFT size

% Pilot Symbols
PS = 1;
% Data Symbols
DS = n_fft - PS;

% OFDM frame length
OFDM = n_fft + n_cp;

% Signal to Noise Ratios for AWGN
snr = [0:1:30];

% Channel modelled by random complex number
h = [randn()+randn()*1i];

% Initialise variables
errors = zeros(size(snr));
ber = zeros(size(snr));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                          TRANSMITTER                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate binary data
binary_data = round(rand(data,1));      

%% Run BER test for 4QAM, 16QAM and 64 QAM
for mod_type = 1:3
    % Quadrature Amplitude Modulation
    if mod_type == 1                                % Modulation Type: 4QAM
        symbols = 2;                                % Bits per symbol  
        pilot = [0 1]';
        % Padding for symbol mapping
        while floor(length(binary_data)/symbols) ~= length(binary_data)/symbols
            binary_data = [binary_data; zeros(1,1)];                                    
        end
        % Padding for subcarriers mapping
        while floor(length(binary_data)/symbols/DS) ~= length(binary_data)/symbols/DS
            binary_data = [binary_data; zeros(1,1)];                                    
        end        
    elseif mod_type == 2                            % Modulation Type: 16QAM
        symbols = 4;
        pilot = [0 0 0 1]';
        % Padding for symbol mapping
        while floor(length(binary_data)/symbols) ~= length(binary_data)/symbols
            binary_data = [binary_data; zeros(1,1)];                                    
        end
        % Padding for subcarriers mapping
        while floor(length(binary_data)/symbols/DS) ~= length(binary_data)/symbols/DS
            binary_data = [binary_data; zeros(1,1)];                                    
        end        
    elseif mod_type == 3                            % Modulation Type: 64QAM
        symbols = 6;
        pilot = [0 0 0 0 0 1]';
        % Padding for symbol mapping
        while floor(length(binary_data)/symbols) ~= length(binary_data)/symbols
            binary_data = [binary_data; zeros(1,1)];                                    
        end
        % Padding for subcarriers mapping
        while floor(length(binary_data)/symbols/DS) ~= length(binary_data)/symbols/DS
            binary_data = [binary_data; zeros(1,1)];                                    
        end
    end
    % Mapping binary data onto constellation maps
    mod_method = 2^symbols;                        
    mod_data = qammod(binary_data,mod_method,'unitaveragepower',true,'inputtype','bit');
    pilot = qammod(pilot,mod_method,'unitaveragepower',true,'inputtype','bit');

    Tx1 = mod_data;
    
    %% Serial to Parallel Conversion
    % serial stream into subcarriers
    Xk1 = reshape(Tx1,DS,length(Tx1)/DS);
    
    %% Pilot Insertion
    Xk1_p = [ones(1,length(Xk1))*pilot;Xk1];
   
    %% Inverse Fast Fourier Transform
    Xn1 = ifft(Xk1_p);

    %% Cyclic prefix
    Xn1_cp = [Xn1((end - n_cp + 1):end,:);Xn1];  

    %% Parallel to Serial Conversion
    xn1 = Xn1_cp(:);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              Channel                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% BER test for different SNR
    for k = 1:length(snr)
        %% Data passing through Channel
        % y time, receiver 
        yn11 = xn1*h(1,1);              % Tx1 to Rx1               

        %% add noise
        dB = snr(k);
        yn11 = awgn(yn11,dB,'measured');

        %% offset added to between transmitter and receiver
        % CP length is 16, input time offset of 15 to match CP length of 16 
        % offset greater than CP length, 15, introduces ISI and increased BER
        % time offset is equal to duration of one symbol length
        delay11 = 0;          % offset between Tx1 and Rx1

        % next OFDM symbol becomes superposition of delayed previous symbol
        % tail of previous symbol added to front of next symbol - ISI
        %% inter-symbol interference
        for i = 0:(length(yn11)/(OFDM))-2
            % ISI between transmitter 1 and receiver 1
            if delay11 ~= 0
                if delay11 == 1
                    yn11(1+OFDM+OFDM*i) = yn11(1+OFDM+OFDM*i) + yn11(OFDM-delay11+OFDM*i);
                else
                    yn11(1+OFDM+OFDM*i:1+OFDM+delay11+OFDM*i) = yn11(1+OFDM+OFDM*i:1+OFDM+delay11+OFDM*i) ...
                    + yn11(OFDM-delay11+OFDM*i:OFDM+OFDM*i);
                end
            end
        end    
        % Time offset transmitter 1 and receiver 1
        if delay11 ~= 0
            if delay11 == 1
                 yn11(1:OFDM) = [zeros(delay11,1);yn11(1:OFDM-delay11)];
            else
                yn11(1:OFDM) = [zeros(delay11+1,1);yn11(1:OFDM-delay11-1)];
            end
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              RECEIVER                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% Serial to Parallel Conversion
        yn1_sp = reshape(yn11,OFDM,length(Xn1_cp));
        
        %% Remove cyclic prefix
        yn1_rcp = yn1_sp((n_cp + 1):end,:);

        %% Discrete fourier transform
        Yk1_block = fft(yn1_rcp);
        %% Channel Estimation
        % because we send pilot symbols, we can estimate the channel. symbols that
        % are known beforehand. So its the received signal divided by the pilot
        % symbol.
        Yk1_DS = Yk1_block(2:64,:);
        Yk1_pilot = Yk1_block(1,:);
        H_hat = Yk1_pilot./pilot;           
        
        %% Channel equalisation
        % zero forcing / Least Squares method
        H_hate = mean(H_hat);                       % average of estimated channel

        Yk1_block2 = Yk1_DS./H_hate;                % cancelling out the channel effects
        
        %% Parallel to serial conversion
        Yk1 = Yk1_block2(:);              

        %% QAM demodulation
        X_demod = qamdemod(Yk1,mod_method,'unitaveragepower',true,'outputtype','bit');
        output = X_demod(:).';
        
        %% Calculating BER
        errors(mod_type,k) = 0;
        for i = 1:length(binary_data)-n_fft*symbols
            if output(i+n_fft*symbols) ~= binary_data(i+n_fft*symbols)
                errors(mod_type,k) = errors(mod_type,k) + 1;
            end
        end
        ber(mod_type,k) = errors(mod_type,k)/length(binary_data);
    end
    
    % Plotting BER vs Eb/No
    semilogy(snr-10*log10(symbols),ber(mod_type,:),'-');
    title('Single transmitter Single Receiver'); legend('4QAM','16QAM','64QAM');
    xlabel('E_b/N_o (dB)'); ylabel('BER'); grid on; hold on
end
