% OFDM two transmitters one receiver transmission
% Name: Benjamin Pham
% Student ID: 25957066
% 
% Task: to simulate a transmission of data using OFDM technique
% blocks to so simulate, modulation block, IFFT, P/S,
% FFT, demodulate, recover signal. then add (2) Noise,(3) Multipath Channel

% Now integretating a second transmitter, before using space time block
% coding, will use IFFT where half the data is zeros so when superimposed
% it will become one full transmission
% the second transmitter will send data using an a range of allocated
% subcarriers
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
% Rb = 100*10^6;y
% B = 20*10^6;

% fft/ifft size
n_fft = 256;                                            % increased fft size to from 64
%pilot length
n_p = 4;

% cyclic prefix size
n_cp = 16;

% snr
snr = [0:1:40];

% attenuation
for a = 1:1                    % change 64 to 1, to produce just a single complex number, equalising works
chngain(a) = rand() + i*rand();
end

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
t_data1 = round(rand(bits,1));                                               % generate data
t_data1 = dec2bin(t_data1);                                                  % changes vector to char type

t_data2 = round(rand(bits,1));                                              
t_data2 = dec2bin(t_data2);                                                 


%% symbol mapping
% 4QAM, 16QAM

reshape_data1 = reshape(t_data1,length(t_data1)/symbols, symbols );                       % reshape into pairings of 2 bits per symbol
reshape_data2 = reshape(t_data2,length(t_data2)/symbols, symbols ); 

decimal_data1 = bin2dec(reshape_data1);                                                   % must be a char vector
decimal_data2 = bin2dec(reshape_data2);      

if mod_value == 1
    mod_data1 = qammod(decimal_data1,4,'unitaveragepower',true);
    mod_data2 = qammod(decimal_data2,4,'unitaveragepower',true);
else 
    mod_data1 = qammod(decimal_data1,16,'unitaveragepower',true);
    mod_data2 = qammod(decimal_data2,16,'unitaveragepower',true);
end 
%figure()    
%scatterplot(mod_data)

%% reshape

X1 = mod_data1;
X2 = mod_data2;

X_blocks1 = reshape(X1,n_fft/2,length(X1)/(n_fft/2));                        % reshape into 128 block
X_blocks2 = reshape(X2,n_fft/2,length(X2)/(n_fft/2));                       

X_blocks_pad1 = [ X_blocks1 ; zeros(size(X_blocks1)) ] ; 
X_blocks_pad2 = [ zeros(size(X_blocks2)) ; X_blocks2 ] ;                                                             

X_blocks = X_blocks_pad1 + X_blocks_pad2;                                    % superimpose the two data blocks from different subcarriers

%% insert pilot symbols 
% pilot symbol is inserted on each subcarrier (block type)
pilots = ones(1,1);

if mod_value == 1
    mod_pilots = qammod(pilots,4,'unitaveragepower',true);
else 
    mod_pilots = qammod(pilots,16,'unitaveragepower',true);
end 
[c,d] = size(X_blocks);

X_block = zeros(c,d+1);
for i2 = 1:n_fft
X_block(i2,:) = [mod_pilots(1,1), X_blocks(i2,:)];  
end 

% mod_pilots(1,2) X_block(i2,(length(X_block)/4+1):length(X_block)/2) ...
% mod_pilots(1,3) X_block(i2,(length(X_block)/2+1):3*(length(X_block)/4)) mod_pilots(1,4) X_block(i2,(3*(length(X_block)/4)+1):end)];

%% IFFT
% Moves the data from the freq domain to the time domain
x = ifft(X_block);                                                          % Inverse fast fourier transform

%% add CP
x_cp = [x(:,(end - n_cp + 1) :end),x];                                      % add CP to end of data                                 
%% Parallel data to Serial stream

x_s = x_cp(:);

%%                          CHANNEL
%% Multipath Channel
% delays and attentuation that affects signal      
H_x = conv(x_s,chngain,'same');
% H_x = x_s;

%% AWGN Noise
for i = 1:length(snr)
    H_noise = awgn(H_x,snr(i),'measured');
    %H_noise = H_x;
%%                          RECEIVER
%% Serial to Parallel
y_p = reshape(H_noise, n_fft, length(H_noise + n_cp - 1)/n_fft);

% remove cp
y_p_rcp = y_p(:,(n_cp + 1):end);

%% FFT
% converts signal from time domain to frequency domain
Y_blocks = fft(y_p_rcp);
Y = reshape(Y_blocks, (length(X1)+length(X2)+n_fft*length(pilots)),1);

%% Channel Estimation
% because we send pilot symbols, we can estimate the channel. symbols that
% are known beforehand. So its the received signal divided by the pilot
% symbol.

Y_hat = zeros(c,d+1);


for i2 = 1:n_fft/2
H_hat1(i2) = Y_blocks(i2,1)./X_block(i2,1);                            % estimating channel using pilot symbols
H_hat2(i2) = Y_blocks(i2+n_fft/2,1)./X_block(i2+n_fft/2,1);  

Y_hat(i2,:) = H_hat1(i2) .* X_block(i2,:);
Y_hat(i2+n_fft/2,:) = H_hat2(i2) .* X_block(i2+n_fft/2,:);     
end
             
% Y_hat1 = H_hat1 .* X_blocks1(1:n_p);                  
% Y_hat2 = H_hat2 .* X_blocks2(length(Y)/2:length(Y)/2+ n_p); 

% % mean squared error
% for m = 1:length(H_hat)
%     se(m) = (abs(Y(m) - Y_hat(m)))^2;
%     mse(m) = sum(se)/length(se);
% end
% mmse(mod_value,i) = min(mse);

 %% equalisation
H_hat_e1 = mean(H_hat1);
H_hat_e2 = mean(H_hat2);

Y_blocks2 = zeros(c,d+1);

Y_blocks2(1:n_fft/2,:) = Y_blocks(1:n_fft/2,:)/H_hat_e1;                     % cancelling out the channel effects?
Y_blocks2(n_fft/2+1:end,:) = Y_blocks(n_fft/2+1:end,:)/H_hat_e2;   

% remove pilot
Y_block = Y_blocks(:,2:end);
Y_blocks2 = Y_blocks2(:,2:end);
%% Demodulate
if mod_value == 1
    y = qamdemod(Y_block,4,'unitaveragepower',true);
else 
    y = qamdemod(Y_block,16,'unitaveragepower',true);
end

received_sym1 = dec2bin(y(1:n_fft/2,:));
received_sym2 = dec2bin(y((n_fft/2+1):end,:));

received_sym = [received_sym1,received_sym2];

received_sig = reshape(received_sym, length(t_data1)+length(t_data2), 1);

errors_t1(mod_value,i) = 0;
for k1 = [1:length(received_sig)/2]
    if received_sig(k1) ~= t_data1(k1)
        errors_t1(mod_value,i) = errors_t1(mod_value,i) + 1;
    end
end

errors_t2(mod_value,i) = 0;
for k2 = [length(received_sig)/2+1:length(received_sig)]
    if received_sig(k2) ~= t_data2(k2-length(received_sig)/2)
        errors_t2(mod_value,i) = errors_t2(mod_value,i) + 1;
    end
end

ber_t1(mod_value,i) = errors_t1(mod_value,i)/length(t_data1);
ber_t2(mod_value,i) = errors_t2(mod_value,i)/length(t_data2);
ber(mod_value,i) = (errors_t1(mod_value,i)+ errors_t2(mod_value,i))/(length(t_data1)+length(t_data2));
%% Demodulate equalistation part
if mod_value == 1
    y2 = qamdemod(Y_blocks2,4,'unitaveragepower',true);
else 
    y2 = qamdemod(Y_blocks2,16,'unitaveragepower',true);
end

received_sym1_eq = dec2bin(y2(1:n_fft/2,:));
received_sym2_eq = dec2bin(y2((n_fft/2+1):end,:));

received_sym_eq = [received_sym1_eq,received_sym2_eq];

received_sig_eq = reshape(received_sym_eq, length(t_data1)+length(t_data2), 1);

errors2_t1(mod_value,i) = 0;

for m1 = [1:length(received_sig_eq)/2]
    if received_sig_eq(m1) ~= t_data1(m1)
        errors2_t1(mod_value,i) = errors2_t1(mod_value,i) + 1;
    end
end
errors2_t2(mod_value,i) = 0;
for m2 = [1:length(received_sig_eq)/2]
    if received_sig_eq(m2+(length(received_sig_eq)/2)) ~= t_data2(m2)
        errors2_t2(mod_value,i) = errors2_t2(mod_value,i) + 1;
    end
end

ber2_t1(mod_value,i) = errors2_t1(mod_value,i)/length(t_data1);
ber2_t2(mod_value,i) = errors2_t2(mod_value,i)/length(t_data2);
ber2(mod_value,i) = (errors2_t1(mod_value,i)+ errors2_t2(mod_value,i))/(length(t_data1)+length(t_data2));
end
%figure()
%semilogy(snr,ber(mod_value,:),'x-',snr,ber2(mod_value,:),'o-');
semilogy(snr,ber2(mod_value,:),'o-');
title('BER vs SNR'); legend('4QAM','16QAM');
xlabel('SNR'); ylabel('BER'); grid on; hold on


end


