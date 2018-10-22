clc; clear all; close all

data = 2^16;                                    % data points

n_fft = 64;                                     % fft size
    
n_cp = 16;                                      % cyclic prefix

snr = [0:1:35];
errors = zeros(size(snr));
ber = zeros(size(snr));
OFDM = n_fft+n_cp;
%%
binary_data = round(rand(data,1));              % generate data

h = [randn()+randn()*1i; randn()+randn()*1i];
%%
for mod_value = 1:3
% bits per symbol
if mod_value == 1
    symbols = 2;   
elseif mod_value == 2
    symbols = 4;
else 
    symbols = 6;
    while floor(length(binary_data)/symbols) ~= length(binary_data)/symbols
    binary_data = [binary_data; zeros(1,1)];                                    %padding for symbol mapping
    end
    while floor(length(binary_data)/symbols/n_fft) ~= length(binary_data)/symbols/n_fft
    binary_data = [binary_data; zeros(6,1)];                                    %padding for subcarriers mapping
    end
end

mod_method = 2^symbols;       
                                                
mod_data = qammod(binary_data,mod_method,'unitaveragepower',true,'inputtype','bit');
%mod_data = mod_data./abs(mod_data);
%% STBC
% [x1 -x2* ; x2 x1*]
STBC = zeros(2*length(mod_data),1);
for i = 0:(length(mod_data)/2-1)
    STBC(4*i+1) = mod_data(2*i+1);
    STBC(4*i+2) = mod_data(2*i+2);
    STBC(4*i+3) = -conj(mod_data(2*i+2));
    STBC(4*i+4) = conj(mod_data(2*i+1));
end

STBC = reshape(STBC.',2,length(mod_data));
%% splitting up the data between the transmitters
Tx1 = STBC(1,:);
Tx2 = STBC(2,:);

% while floor(length(Tx1)/64) ~= length(Tx1)/64
%     Tx1 = [Tx1 zeros(1,1)];
%     Tx2 = [Tx2 zeros(1,1)];
% end


%% IFFT
Xk1 = reshape(Tx1,n_fft,length(Tx1)/n_fft);
Xk2 = reshape(Tx2,n_fft,length(Tx2)/n_fft);

Xn1 = ifft(Xk1);
Xn2 = ifft(Xk2);
%% Cyclic prefix
Xn1_cp = [Xn1((end - n_cp + 1):end,:);Xn1];  
Xn2_cp = [Xn2((end - n_cp + 1):end,:);Xn2]; 

%% P/S
xn1 = Xn1_cp(:);
xn2 = Xn2_cp(:);
%% Channel

for k = 1:length(snr)
yn11 = xn1*h(1);            % y time, receiver 
yn21 = xn2*h(2);

% add noise

yn11 = awgn(yn11,snr(k),'measured');
yn21 = awgn(yn21,snr(k),'measured');

%% Delay added
delay11 = 19;
delay21 = 16;
% adding delay to first transmitter


% yn11 = [zeros(delay,1);yn11];
% yn21 = [yn21; zeros(delay,1)];
for i = 0:(length(yn11)/(OFDM))-2
    
    if delay11 ~= 0
    yn11(1+OFDM+OFDM*i:1+OFDM+delay11+OFDM*i) = yn11(1+OFDM+OFDM*i:1+OFDM+delay11+OFDM*i) ...
        + yn11(OFDM-delay11+OFDM*i:OFDM+OFDM*i);
    end
    if delay21 ~= 0
    yn21(1+OFDM+OFDM*i:1+OFDM+delay21+OFDM*i) = yn21(1+OFDM+OFDM*i:1+OFDM+delay21+OFDM*i) ...
        + yn21(OFDM-delay21+OFDM*i:OFDM+OFDM*i);
    end

end

if delay11 ~= 0
yn11(1:OFDM) = [zeros(delay11+1,1);yn11(1:OFDM-delay11-1)];      
end
if delay21 ~= 0
yn21(1:OFDM) = [zeros(delay21+1,1);yn21(1:OFDM-delay21-1)];  
end



% superposition of two received signals at receiver
yn1 = yn11 + yn21;                    % y(1) = h1*x1 + h2*x2 ; y(2) = -h2* x1* + h1 * x2*
% yn1 = yn1(1:end-delay);
%% S/P
yn1_sp = reshape(yn1,n_fft+n_cp,length(Xn1_cp));

%% remove cyclic prefix
yn1_rcp = yn1_sp((n_cp + 1):end,:);


%% DFT
Yk1_block = fft(yn1_rcp);


Yk1 = Yk1_block(:);              % transform 64 subchannels back into 1 stream


%% STBC decoding
abs_h = sum(sum(abs(h).^2));
% assume perfect channel estimation
H = [ conj(h(1)) , h(2),; ...         % pseudo inverse
    conj(h(2)), -h(1)]./ abs_h;

X_hat1 = zeros(length(Tx1)/2,1);
X_hat2 = zeros(length(Tx2)/2,1);

X = zeros(length(Tx1)/2,1);

for i = 0:length(Yk1)/2-1
Yk1(2*i+2) = conj(Yk1(2*i+2));

X_hat1(i+1) = H(1,1)*Yk1(2*i+1) + H(1,2)*Yk1(2*i+2);       
X_hat2(i+1) = H(2,1)*Yk1(2*i+1) + H(2,2)*Yk1(2*i+2);

X(2*i+1) = X_hat1(i+1);
X(2*i+2) = X_hat2(i+1);

end

%% demodulate
X_demod = qamdemod(X,mod_method,'unitaveragepower',true,'outputtype','bit');
%output = reshape(X_demod.',size(binary_data));
output = X_demod(:).';
%%

errors(mod_value,k) = 0;
for i = 1:length(binary_data)-n_fft*symbols
    if X_demod(i+n_fft*symbols) ~= binary_data(i+n_fft*symbols)
        errors(mod_value,k) = errors(mod_value,k) + 1;
    end
end
ber(mod_value,k) = errors(mod_value,k)/length(binary_data);

end

semilogy(snr-10*log10(symbols),ber(mod_value,:),'-');
title('STBC for 2x1 system'); legend('4QAM','16QAM','64QAM');
xlabel('SNR'); ylabel('BER'); grid on; hold on
end