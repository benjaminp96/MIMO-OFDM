%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Student:      Benjamin Pham
%   ID:           25957066
%   
%   Version:        Final
%   Changes:        
%   Description:    a two user one receiver OFDM System
%                                   
%   Simulates a multi-user OFDM system. This system includes the
%   Transmitter: 
%   QAM Modulation, S/P converter, IFFT, Add Cyclic Prefix, 
%   Sub-carrier allocation, P/S converter
%   Channel: 
%   Frequency flat fading Channel modelled by Complex number, AWGN, 
%   Time Offset
%   Receiver: 
%   S/P converter, Sub-carrier Demux, Removal of Cyclic Prefix, DFT, 
%   P/S converter, QAM Demodulation
%
%   Simulation performs a BER test vs SNR when various time offsets are
%   inputted by sending binary data through the system and measuring the
%   amount of errors at the receiver. Time offsets can be changed beginning
%   on line 140 for different channels.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all;close all; warning off
%% Input parameters
% data points
data = 2^16;                                   

% fft size
n_fft = 64;                                     

% cyclic prefix length
n_cp = 0.25*n_fft;                                  % CP length 25% of FFT size

% OFDM frame length
OFDM = 2*n_fft+2*n_cp;

% Signal to Noise Ratios for AWGN
snr = [0:1:30];

% Channel modelled by random complex number
% 1 for same channels or any other number different channels
channel = 1;
% Simulate users with same channel or different channels
if channel == 1
    h = [randn() + randn()*1i];                                  
    h = [h h];
else
    h = [randn()+randn()*1i , randn()+randn()*1i]; 
end

% Initialise variables
errors1 = zeros(size(snr));
ber1 = zeros(size(snr));
errors2 = zeros(size(snr));
ber2 = zeros(size(snr));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                          TRANSMITTER                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate binary data
binary_data1 = round(rand(data,1));             % User 1
binary_data2 = round(rand(data,1));             % User 2

%% Run BER test for 4QAM, 16QAM and 64 QAM
for mod_type = 1:3
    % Quadrature Amplitude Modulation
    if mod_type == 1                                % Modulation Type: 4QAM
        symbols = 2;                                % Bits per symbol  
        % Padding for symbol mapping
        while floor(length(binary_data1)/symbols) ~= length(binary_data1)/symbols
        binary_data1 = [binary_data1; zeros(1,1)];                                    
        end
        while floor(length(binary_data1)/symbols/(2*n_fft)) ~= length(binary_data1)/symbols/(2*n_fft)
        binary_data1 = [binary_data1; zeros(1,1)];                                    
        end
        % Padding for subcarriers mapping
        while floor(length(binary_data2)/symbols) ~= length(binary_data2)/symbols
        binary_data2 = [binary_data2; zeros(1,1)];                                    
        end
        while floor(length(binary_data2)/symbols/(2*n_fft)) ~= length(binary_data2)/symbols/(2*n_fft)
        binary_data2 = [binary_data2; zeros(1,1)];                                    
        end       
    elseif mod_type == 2                            % Modulation Type: 16QAM
        symbols = 4;
        % Padding for symbol mapping
        while floor(length(binary_data1)/symbols) ~= length(binary_data1)/symbols
        binary_data1 = [binary_data1; zeros(1,1)];                                    
        end
        while floor(length(binary_data1)/symbols/(2*n_fft)) ~= length(binary_data1)/symbols/(2*n_fft)
        binary_data1 = [binary_data1; zeros(1,1)];                                    
        end
        % Padding for subcarriers mapping
        while floor(length(binary_data2)/symbols) ~= length(binary_data2)/symbols
        binary_data2 = [binary_data2; zeros(1,1)];                                    
        end
        while floor(length(binary_data2)/symbols/(2*n_fft)) ~= length(binary_data2)/symbols/(2*n_fft)
        binary_data2 = [binary_data2; zeros(1,1)];                                    
        end     
    elseif mod_type == 3                            % Modulation Type: 64QAM
        symbols = 6;
        % Padding for symbol mapping
        while floor(length(binary_data1)/symbols) ~= length(binary_data1)/symbols
        binary_data1 = [binary_data1; zeros(1,1)];                                    
        end
        while floor(length(binary_data1)/symbols/(2*n_fft)) ~= length(binary_data1)/symbols/(2*n_fft)
        binary_data1 = [binary_data1; zeros(1,1)];                                    
        end
        % Padding for subcarriers mapping
        while floor(length(binary_data2)/symbols) ~= length(binary_data2)/symbols
        binary_data2 = [binary_data2; zeros(1,1)];                                    
        end
        while floor(length(binary_data2)/symbols/(2*n_fft)) ~= length(binary_data2)/symbols/(2*n_fft)
        binary_data2 = [binary_data2; zeros(1,1)];                                    
        end
    end
    % Mapping binary data onto constellation maps
    mod_method = 2^symbols;       
    mod_data1 = qammod(binary_data1,mod_method,'unitaveragepower',true,'inputtype','bit');
    mod_data2 = qammod(binary_data2,mod_method,'unitaveragepower',true,'inputtype','bit');

    %% splitting data up between the users
    Tx1 = mod_data1;
    Tx2 = mod_data2;

    %% Serial to Parallel Conversion
    % serial stream into 64 subcarriers
    Xk1 = reshape(Tx1,n_fft,length(Tx1)/(n_fft)); % 64 sub carriers, per user
    Xk2 = reshape(Tx2,n_fft,length(Tx2)/(n_fft));

    %% Inverse Fast Fourier Transform
    Xn1 = ifft(Xk1);
    Xn2 = ifft(Xk2);
    
    %% Cyclic prefix
    Xn1_cp = [Xn1((end - n_cp + 1):end,:);Xn1];  
    Xn2_cp = [Xn2((end - n_cp + 1):end,:);Xn2];  
  
    %% Sub-carrier allocation    
    Xn1_pad = [ Xn1_cp ; zeros(size(Xn1_cp)) ];
    Xn2_pad = [ zeros( size(Xn2_cp)) ; Xn2_cp ];

    %% Parallel to Serial Conversion
    xn1 = Xn1_pad(:);
    xn2 = Xn2_pad(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              Channel                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BER test for different SNR
    for k = 1:length(snr)
    %% Data passing through Channel
        yn11 = xn1*h(1);           
        yn21 = xn2*h(2);
        
        %% add noise
        db = snr(k);
        yn11 = awgn(yn11,db,'measured');
        yn21 = awgn(yn21,db,'measured');
        
        %% offset added to between User and receiver
        % CP length is 16, input time offset of 15 to match CP length of 16 
        % offset greater than CP length, 15, introduces ISI and increased BER
        % time offset is equal to duration of one symbol length
        delay11 = 0;          % offset between User1 and Rx1
        delay21 = 0;          % offset between User2 and Rx1

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
            % ISI between transmitter 1 and receiver 2            
            if delay21 ~= 0
                if delay21 == 1
                    yn21(1+OFDM+OFDM*i) = yn21(1+OFDM+OFDM*i) + yn21(OFDM-delay21+OFDM*i);
                else                
                yn21(1+OFDM+OFDM*i:1+OFDM+delay21+OFDM*i) = yn21(1+OFDM+OFDM*i:1+OFDM+delay21+OFDM*i) ...
                    + yn21(OFDM-delay21+OFDM*i:OFDM+OFDM*i);
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
        % Time offset transmitter 1 and receiver 2        
        if delay21 ~= 0
            if delay21 == 1
                 yn21(1:OFDM) = [zeros(delay21,1);yn21(1:OFDM-delay21)];
            else
                yn21(1:OFDM) = [zeros(delay21+1,1);yn21(1:OFDM-delay21-1)];
            end
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              RECEIVER                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        % superposition of two received signals at receiver
        %% y(1) = h1*x1 + h2*x2
        yn = yn11 + yn21;               % combined signal at receiver 1  
        
        %% Serial to Parallel Conversion
        yn_sub = reshape(yn,OFDM,size(Xn1,2));
        
        %% Sub-Carrier Demux
        Yn_sub1 = yn_sub(1:80,:);                              
        Yn_sub2 = yn_sub(81:160,:);

        %% Remove cyclic prefix
        yn_sub_rcp1 = Yn_sub1((n_cp+1):end,:);
        yn_sub_rcp2 = Yn_sub2((n_cp+1):end,:);


        %% Discrete fourier transform
        Yk_block1 = fft(yn_sub_rcp1);
        Yk_block2 = fft(yn_sub_rcp2);
        
        %% Serial to parallel conversion
        Yk1_s = Yk_block1(:);             
        Yk2_s = Yk_block2(:);       

        %% Channel estimation
        % assume perfect channel estimation through pilot symbols and using
        % zero forcing
        Xhat1 = Yk1_s./h(1);
        Xhat2 = Yk2_s./h(2);

        %% QAM demodulation
        X_demod1 = qamdemod(Xhat1,mod_method,'unitaveragepower',true,'outputtype','bit');
        X_demod2 = qamdemod(Xhat2,mod_method,'unitaveragepower',true,'outputtype','bit');

        output1 = X_demod1(:).';
        output2 = X_demod2(:).';
        
        %% Calculating BER
        % BER for User 1
        errors1(mod_type,k) = 0;
        for i = 1:length(binary_data1)-n_fft*symbols
            if output1(i+n_fft*symbols) ~= binary_data1(i+n_fft*symbols)
                errors1(mod_type,k) = errors1(mod_type,k) + 1;
            end
        end
        ber1(mod_type,k) = errors1(mod_type,k)/length(binary_data1);
        
        % BER for User 2
        errors2(mod_type,k) = 0;
        for i = 1:length(binary_data2)-n_fft*symbols
            if output2(i+n_fft*symbols) ~= binary_data2(i+n_fft*symbols)
                errors2(mod_type,k) = errors2(mod_type,k) + 1;
            end
        end
        ber2(mod_type,k) = errors2(mod_type,k)/length(binary_data2);

    end
    
    % Plotting BER vs Eb/No
    semilogy(snr-10*log10(symbols),ber1(mod_type,:),'-',snr-10*log10(symbols),ber2(mod_type,:),'-')
    title('Two users One receiver'); legend('4QAM user1','4QAM user2','16QAM user1','16QAM user2','64QAM user1','64QAM user2');
    xlabel('E_b/N_o (dB)'); ylabel('BER'); grid on;hold on
end
