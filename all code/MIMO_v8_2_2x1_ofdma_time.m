clc; clear all; close all; warning off

data = 2^16;                                    % data points

n_fft = 128;                                     % fft size
    
n_cp = 16;                                      % cyclic prefix

snr = [0:1:35];

errors1 = zeros(size(snr));
ber1 = zeros(size(snr));
errors2 = zeros(size(snr));
ber2 = zeros(size(snr));
%%
binary_data1 = round(rand(data,1));              % generate data
binary_data2 = round(rand(data,1));              % generate data


% h = [randn()+randn()*1i; randn()+randn()*1i];

h = [randn()+randn()*1i];
h = [h h];

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
Xk1 = reshape(Tx1,n_fft/2,length(Tx1)/(n_fft/2)); % 64 sub carriers, 
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
yn11 = Xn1_cp*h(1);            % y time, receiver 
yn21 = Xn2_cp*h(2);

% add noise

yn11 = awgn(yn11,snr(k),'measured');
yn21 = awgn(yn21,snr(k),'measured');

% superposition of two received signals at receiver
yn = yn11 + yn21;                    % y(1) = h1*x1 + h2*x2 ; y(2) = -h2* x1* + h1 * x2*


%% remove cyclic prefix
yn_rcp = yn(:,(n_cp + 1):end);


%% DFT
Yk_block = fft(yn_rcp);

% decode sub bands to seperate users
Yk1 = Yk_block(1:64,:);                               % higher order decision block
Yk2 = Yk_block(65:128,:);

Yk1_s = reshape(Yk1, 1, length(Tx1));              % transform 64 subchannels back into 1 stream
Yk2_s = reshape(Yk2, 1, length(Tx2));            



%% Channel estimation
% assume perfect channel estimation, so we know what h1 and h2 are

Xhat1 = Yk1_s./h(1);
Xhat2 = Yk2_s./h(2);

%% demodulate
X_demod1 = qamdemod(Xhat1,mod_method,'unitaveragepower',true,'outputtype','bit');
X_demod2 = qamdemod(Xhat2,mod_method,'unitaveragepower',true,'outputtype','bit');

output1 = X_demod1(:).';
output2 = X_demod2(:).';
%%

errors1(mod_value,k) = 0;
for i = 1:length(binary_data1)
    if output1(i) ~= binary_data1(i)
        errors1(mod_value,k) = errors1(mod_value,k) + 1;
    end
end
ber1(mod_value,k) = errors1(mod_value,k)/length(binary_data1);


errors2(mod_value,k) = 0;
for i = 1:length(binary_data2)
    if output2(i) ~= binary_data2(i)
        errors2(mod_value,k) = errors2(mod_value,k) + 1;
    end
end
ber2(mod_value,k) = errors2(mod_value,k)/length(binary_data2);


end
semilogy(snr,ber1(mod_value,:),'x-',snr,ber2(mod_value,:),'o-')
%semilogy(snr,ber2(mod_value,:),'o-');
title('OFDMA for 2x1 system'); legend('4QAM user1','4QAM user2','16QAM user1','16QAM user2','64QAM user1','64QAM user2');
xlabel('SNR'); ylabel('BER'); grid on;hold on

end