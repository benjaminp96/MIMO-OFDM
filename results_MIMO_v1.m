%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Student:      Benjamin Pham
%   ID:           25957066
%   
%   Version:        Final
%   Changes:        
%   Description:    a 2x2 MIMO OFDM System
%                                   
%   Simulates a 2x2 MIMO OFDM system. This system includes the
%   Transmitter: 
%   QAM Modulation, S/P converter, Space-Time Block encoding
%   IFFT, Add Cyclic Prefix, P/S converter
%   Channel: 
%   Frequency flat fading Channel modelled by Complex number, AWGN, 
%   Time Offset
%   Receiver: 
%   S/P converter, Removal of Cyclic Prefix, DFT, P/S converter
%   Space-time block decoding, QAM Demodulation
%
%   Simulation performs a BER test vs SNR when various time offsets are
%   inputted by sending binary data through the system and measuring the
%   amount of errors at the receiver. Time offsets can be changed beginning
%   on line 130 for different channels.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all; warning off

%% Input parameters
% data points
data = 2^16;                                   

% fft size
n_fft = 64;                                     

% cyclic prefix length
n_cp = 0.25*n_fft;                                  % CP length 25% of FFT size

% OFDM frame length
OFDM = n_fft + n_cp;

% Signal to Noise Ratios for AWGN
snr = [0:1:30];

% Channel modelled by random complex number
h = [randn()+randn()*1i, randn()+randn()*1i; randn()+randn()*1i, randn()+randn()*1i];

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
        % Padding for symbol mapping
        while floor(length(binary_data)/symbols) ~= length(binary_data)/symbols
            binary_data = [binary_data; zeros(1,1)];                                    
        end
        % Padding for subcarriers mapping
        while floor(length(binary_data)/symbols/n_fft) ~= length(binary_data)/symbols/n_fft
            binary_data = [binary_data; zeros(1,1)];                                    
        end        
    elseif mod_type == 2                            % Modulation Type: 16QAM
        symbols = 4;
        % Padding for symbol mapping
        while floor(length(binary_data)/symbols) ~= length(binary_data)/symbols
            binary_data = [binary_data; zeros(1,1)];                                    
        end
        % Padding for subcarriers mapping
        while floor(length(binary_data)/symbols/n_fft) ~= length(binary_data)/symbols/n_fft
            binary_data = [binary_data; zeros(1,1)];                                    
        end        
    elseif mod_type == 3                            % Modulation Type: 64QAM
        symbols = 6;
        % Padding for symbol mapping
        while floor(length(binary_data)/symbols) ~= length(binary_data)/symbols
            binary_data = [binary_data; zeros(1,1)];                                    
        end
        % Padding for subcarriers mapping
        while floor(length(binary_data)/symbols/n_fft) ~= length(binary_data)/symbols/n_fft
            binary_data = [binary_data; zeros(1,1)];                                    
        end
    end
    % Mapping binary data onto constellation maps
    mod_method = 2^symbols;                        
    mod_data = qammod(binary_data,mod_method,'unitaveragepower',true,'inputtype','bit');

    %% Space-Time Block Encoding
    % intialiase matrix
    STBC = zeros(2*length(mod_data),1);
    % [x1 -x2* ; x2 x1*] alogritm 
    % first column is first time slot, next column is next time slot
    for i = 0:(length(mod_data)/2-1)
        STBC(4*i+1) = mod_data(2*i+1);
        STBC(4*i+2) = mod_data(2*i+2);
        STBC(4*i+3) = -conj(mod_data(2*i+2));
        STBC(4*i+4) = conj(mod_data(2*i+1));
    end
    %% Splitting up the data between the transmitters
    STBC = reshape(STBC.',2,length(mod_data));
    Tx1 = STBC(1,:);                % Data to be sent from Tx1
    Tx2 = STBC(2,:);                % Data to be sent from Tx2

    %% Serial to Parallel Conversion
    % serial stream into 64 subcarriers
    Xk1 = reshape(Tx1,n_fft,length(Tx1)/n_fft);
    Xk2 = reshape(Tx2,n_fft,length(Tx2)/n_fft);
   
    %% Inverse Fast Fourier Transform
    Xn1 = ifft(Xk1);
    Xn2 = ifft(Xk2);

    %% Cyclic prefix
    Xn1_cp = [Xn1((end - n_cp + 1):end,:);Xn1];  
    Xn2_cp = [Xn2((end - n_cp + 1):end,:);Xn2]; 

    %% Parallel to Serial Conversion
    xn1 = Xn1_cp(:);
    xn2 = Xn2_cp(:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              Channel                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% BER test for different SNR
    for k = 1:length(snr)
        %% Data passing through Channel
        % y time, receiver 
        yn11 = xn1*h(1,1);              % Tx1 to Rx1               
        yn21 = xn2*h(2,1);              % Tx2 to Rx1
        yn12 = xn1*h(1,2);              % Tx1 to Rx2
        yn22 = xn2*h(2,2);              % Tx2 to Rx2

        %% add noise
        dB = snr(k);
        yn11 = awgn(yn11,dB,'measured');
        yn21 = awgn(yn21,dB,'measured');
        yn12 = awgn(yn12,dB,'measured');
        yn22 = awgn(yn22,dB,'measured');

        %% offset added to between transmitter and receiver
        % CP length is 16, input time offset of 15 to match CP length of 16 
        % offset greater than CP length, 15, introduces ISI and increased BER
        % time offset is equal to duration of one symbol length
        delay11 = 0;          % offset between Tx1 and Rx1
        delay21 = 0;          % offset between Tx2 and Rx1
        delay12 = 0;          % offset between Tx1 and Rx2
        delay22 = 0;          % offset between Tx2 and Rx2

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
            if delay12 ~= 0
                if delay12 == 1
                    yn12(1+OFDM+OFDM*i) = yn12(1+OFDM+OFDM*i) + yn12(OFDM-delay11+OFDM*i);
                else                
                yn12(1+OFDM+OFDM*i:1+OFDM+delay12+OFDM*i) = yn12(1+OFDM+OFDM*i:1+OFDM+delay12+OFDM*i) ...
                    + yn12(OFDM-delay12+OFDM*i:OFDM+OFDM*i);
                end
            end
            % ISI between transmitter 2 and receiver 1            
            if delay21 ~= 0
                if delay21 == 1
                    yn21(1+OFDM+OFDM*i) = yn21(1+OFDM+OFDM*i) + yn21(OFDM-delay11+OFDM*i);
                else                
                yn21(1+OFDM+OFDM*i:1+OFDM+delay21+OFDM*i) = yn21(1+OFDM+OFDM*i:1+OFDM+delay21+OFDM*i) ...
                    + yn21(OFDM-delay21+OFDM*i:OFDM+OFDM*i);
                end
            end
            % ISI between transmitter 2 and receiver 2            
            if delay22 ~= 0
                if delay22 == 1
                    yn22(1+OFDM+OFDM*i) = yn22(1+OFDM+OFDM*i) + yn22(OFDM-delay11+OFDM*i);
                else                
                yn22(1+OFDM+OFDM*i:1+OFDM+delay22+OFDM*i) = yn22(1+OFDM+OFDM*i:1+OFDM+delay22+OFDM*i) ...
                    + yn22(OFDM-delay22+OFDM*i:OFDM+OFDM*i);
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
        if delay12 ~= 0
            if delay12 == 1
                 yn12(1:OFDM) = [zeros(delay12,1);yn12(1:OFDM-delay12)];
            else
                yn12(1:OFDM) = [zeros(delay12+1,1);yn12(1:OFDM-delay12-1)];
            end
        end
        % Time offset transmitter 2 and receiver 1        
        if delay21 ~= 0
            if delay21 == 1
                 yn21(1:OFDM) = [zeros(delay21,1);yn11(1:OFDM-delay21)];
            else
                yn21(1:OFDM) = [zeros(delay21+1,1);yn11(1:OFDM-delay21-1)];
            end
        end
        % Time offset transmitter 2 and receiver 2        
        if delay22 ~= 0
            if delay22 == 1
                 yn22(1:OFDM) = [zeros(delay22,1);yn22(1:OFDM-delay22)];
            else
                yn22(1:OFDM) = [zeros(delay22+1,1);yn22(1:OFDM-delay22-1)];
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              RECEIVER                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % superposition of two received signals at receiver       
        %% y(1) = h1*x1 + h2*x2 ; y(2) = -h2* x1* + h1 * x2*
        yn1 = yn11 + yn21;              % combined signal at receiver 1                
        yn2 = yn12 + yn22;              % combined signal at receiver 2

        %% Serial to Parallel Conversion
        yn1_sp = reshape(yn1,OFDM,size(Xn1_cp,2));
        yn2_sp = reshape(yn2,OFDM,size(Xn2_cp,2));
        
        %% Remove cyclic prefix
        yn1_rcp = yn1_sp((n_cp + 1):end,:);
        yn2_rcp = yn2_sp((n_cp + 1):end,:);

        %% Discrete fourier transform
        Yk1_block = fft(yn1_rcp);
        Yk2_block = fft(yn2_rcp);
        
        %% Serial to parallel conversion
        Yk1 = Yk1_block(:);              
        Yk2 = Yk2_block(:);

        %% Space-Time Block Decoding
        % Assume perfect channel estimation        
        abs_h = sum(sum(abs(h).^2));
        H = [ conj(h(1,1)) , conj(h(1,2)), h(2,1), h(2,2) ; ...         % pseudo inverse
            conj(h(2,1)), conj(h(2,2)), -h(1,1) , -h(1,2)]./ abs_h;
        
        % initialise matrix
        X_hat1 = zeros(length(Tx1)/2,1);
        X_hat2 = zeros(length(Tx2)/2,1);
        X = zeros(length(Tx1)/2,1);

        for i = 0:length(Yk1)/2-1
            Yk1(2*i+2) = conj(Yk1(2*i+2));
            Yk2(2*i+2) = conj(Yk2(2*i+2));

            X_hat1(i+1) = H(1,1)*Yk1(2*i+1) + H(1,2)*Yk2(2*i+1) + H(1,3)*Yk1(2*i+2) + H(1,4)*Yk2(2*i+2);       
            X_hat2(i+1) = H(2,1)*Yk1(2*i+1) + H(2,2)*Yk2(2*i+1) + H(2,3)*Yk1(2*i+2) + H(2,4)*Yk2(2*i+2);

            X(2*i+1) = X_hat1(i+1);
            X(2*i+2) = X_hat2(i+1);
        end

        %% QAM demodulation
        X_demod = qamdemod(X,mod_method,'unitaveragepower',true,'outputtype','bit');
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
    title('Two transmitter Two Receiver'); legend('4QAM','16QAM','64QAM');
    xlabel('E_b/N_o (dB)'); ylabel('BER'); grid on; hold on
end
