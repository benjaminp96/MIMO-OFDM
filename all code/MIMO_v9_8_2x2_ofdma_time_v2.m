clc; clear all; close all; warning off

data = 2^16;                                    % data points

n_fft = 128;                                     % fft size
    
n_cp = 16;                                      % cyclic prefix

snr = [0:1:35];

% same channel for each stream
%%
binary_data1 = round(rand(data,1));              % generate data
binary_data2 = round(rand(data,1));              % generate data

% h = [randn()+randn()*1i, randn()+randn()*1i; randn()+randn()*1i, randn()+randn()*1i]; 
% h = [randn()+randn()*1i; randn()+randn()*1i];       % using very different channels for the each of the two
% h = [h h];                                             % users, 
h = [randn()+randn()*1i];
h = [h h; h h];


%%
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

%mod_data = mod_data./abs(mod_data);

%% splitting up the data between the transmitters
Tx1 = mod_data1;
Tx2 = mod_data2;


%% subcarriers
Xk1 = reshape(Tx1,n_fft/2,length(Tx1)/(n_fft/2));
Xk2 = reshape(Tx2,n_fft/2,length(Tx2)/(n_fft/2));

Xk1_pad = [ Xk1 ; zeros(size(Xk1)) ];
Xk2_pad = [ zeros( size(Xk2)) ; Xk2 ];


%% ifft
Xn1 = ifft(Xk1_pad);
Xn2 = ifft(Xk2_pad);

%% Cyclic prefix
Xn1_cp = [Xn1(:,(end - n_cp + 1):end),Xn1];  
Xn2_cp = [Xn2(:,(end - n_cp + 1):end),Xn2];  


%% Channel
% h = [randn()+randn()*1i; randn()+randn()*1i];
for k = 1:length(snr)
yn11 = Xn1_cp*h(1,1);            % y time, receiver 
yn21 = Xn2_cp*h(2,1);

yn12 = Xn1_cp*h(1,2);            % y time, receiver 
yn22 = Xn2_cp*h(2,2);

% add noise

yn11 = awgn(yn11,snr(k),'measured');
yn21 = awgn(yn21,snr(k),'measured');

yn12 = awgn(yn12,snr(k),'measured');
yn22 = awgn(yn22,snr(k),'measured');

% superposition of two received signals at receiver
yn1 = yn11 + yn21;                    % y(1) = h1*x1 + h2*x2 ; y(2) = -h2* x1* + h1 * x2*
yn2 = yn12 + yn22; 

%% remove cyclic prefix
yn1_rcp = yn1(:,(n_cp + 1):end);
yn2_rcp = yn2(:,(n_cp + 1):end);


%% DFT
Yk1_block = fft(yn1_rcp);
Yk2_block = fft(yn2_rcp);

% decode sub bands to seperate users
Yk11 = Yk1_block(1:64,:);                               % higher order decision block
Yk21 = Yk1_block(65:128,:);

Yk12 = Yk2_block(1:64,:);                           
Yk22 = Yk2_block(65:128,:);


Yk11_s = reshape(Yk11, 1, length(Tx1));              % transform 64 subchannels back into 1 stream
Yk21_s = reshape(Yk21, 1, length(Tx2));        

Yk12_s = reshape(Yk12, 1, length(Tx1));              % transform 64 subchannels back into 1 stream
Yk22_s = reshape(Yk22, 1, length(Tx2));  



%% Channel estimation
% assume perfect channel estimation, so we know what h1 and h2 are

Xhat11 = Yk11_s./h(1,1);
Xhat21 = Yk21_s./h(2,1);

Xhat12 = Yk12_s./h(1,2);
Xhat22 = Yk22_s./h(2,2);

%% demodulate
X_demod11 = qamdemod(Xhat11,mod_method,'unitaveragepower',true,'outputtype','bit');
X_demod21 = qamdemod(Xhat21,mod_method,'unitaveragepower',true,'outputtype','bit');

output11 = X_demod11(:).';
output21 = X_demod21(:).';

X_demod12 = qamdemod(Xhat11,mod_method,'unitaveragepower',true,'outputtype','bit');
X_demod22 = qamdemod(Xhat21,mod_method,'unitaveragepower',true,'outputtype','bit');

output12 = X_demod12(:).';
output22 = X_demod22(:).';
%%

errors11(mod_value,k) = 0;
for i = 1:length(binary_data1)
    if output11(i) ~= binary_data1(i)
        errors11(mod_value,k) = errors11(mod_value,k) + 1;
    end
end
ber11(mod_value,k) = errors11(mod_value,k)/length(binary_data1);


errors21(mod_value,k) = 0;
for i = 1:length(binary_data2)
    if output21(i) ~= binary_data2(i)
        errors21(mod_value,k) = errors21(mod_value,k) + 1;
    end
end
ber21(mod_value,k) = errors21(mod_value,k)/length(binary_data2);

errors12(mod_value,k) = 0;
for i = 1:length(binary_data1)
    if output12(i) ~= binary_data1(i)
        errors12(mod_value,k) = errors12(mod_value,k) + 1;
    end
end
ber12(mod_value,k) = errors12(mod_value,k)/length(binary_data1);


errors22(mod_value,k) = 0;
for i = 1:length(binary_data2)
    if output22(i) ~= binary_data2(i)
        errors22(mod_value,k) = errors22(mod_value,k) + 1;
    end
end
ber22(mod_value,k) = errors22(mod_value,k)/length(binary_data2);


end

end
semilogy(snr-10*log10(symbols),ber11,'x-',snr-10*log10(symbols),ber21,'^-');
%semilogy(snr,ber11,'x-',snr,ber21,'+-',snr,ber12,'^-',snr,ber22,'v-');
%legend('4QAM user1 Rx1','16QAM user1 Rx1','64QAM user1 Rx1','4QAM user2 Rx1','16QAM user2 Rx1','64QAM user2 Rx1')
%semilogy(snr,ber2(mod_value,:),'o-');
title('OFDMA for 2x2 system receiver 1'); 
legend('4QAM user1','16QAM user1','64QAM user1','4QAM user2','16QAM user2','64QAM user2')
%'4QAM user1 Rx2','16QAM user1 Rx2','64QAM user1 Rx2','4QAM user2 Rx2','16QAM user2 Rx2','64QAM user2 Rx2');
xlabel('SNR'); ylabel('BER'); grid on
figure()
semilogy(snr,ber12,'+-',snr,ber22,'v-')
%legend('4QAM user1 Rx2','16QAM user1 Rx2','64QAM user1 Rx2','4QAM user2 Rx2','16QAM user2 Rx2','64QAM user2 Rx2');
title('OFDMA for 2x2 system receiver 2'); 
legend('4QAM user1','16QAM user1','64QAM user1','4QAM user2','16QAM user2','64QAM user2');
xlabel('SNR'); ylabel('BER'); grid on
