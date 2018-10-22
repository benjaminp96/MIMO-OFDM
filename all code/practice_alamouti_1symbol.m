clc;clear all;close all;
% generate data
N = 4; %1000 bits
symbols = 2;
%n_fft = 64;

data = dec2bin(round(rand(1,N)));
% modulate data using QPSK
data_sym = reshape(data,N/symbols,symbols);     % [ group the bits into bits per symbol]
data_dec = bin2dec(data_sym);                   % convert to decimal
mod_data = qammod(data_dec,4,'unitaveragepower',true);          % 4QAM


% alamouti STBC
% [x1 -x2*; x2 x1*]
STBC = zeros(2,N/symbols);                               
STBC(:,1:2:end) = reshape(mod_data,2,N/(symbols*2));                        % creating [x1 ; x2]
STBC(:,2:2:end) = [-1;1].*flipud(reshape(conj(mod_data),2,N/(symbols*2)));  % creating [-x2* ; x1*]

% channel 
% x + j*y           [h1 ; h2]
h = [randn()+randn()*1i; randn()+rand()*1i];
h(1) = h(1)/norm(h(1));
h(2) = h(2)/norm(h(2));

Tx1 = ifft(STBC(1,:));
Tx2 = ifft(STBC(2,:));

y1 = conv(Tx1,h(1),'same');
y2 = conv(Tx2,h(2),'same');

% y1 = Tx1*h(1);
% y2 = Tx2*h(2);

% superposition of two received signals at receiver
y = y1 + y2;                    % y(1) = h1*x1 + h2*x2 ; y(2) = -h2* x1* + h1 * x2*

% serial to parallel
% y_sub = reshape(y,n_fft,length(STBC)/n_fft);

% DFT 
% Y_blocks  = fft(y_sub);
% lets say we estimated the channel and we know what h1 and h2 are

y = fft(y);

h11 = conj(h(1))/(abs(h(1))^2+abs(h(2))^2);
h21 = h(2)/(abs(h(1))^2+abs(h(2))^2);

h12 = conj(h(2))/(abs(h(1))^2+abs(h(2))^2);
h22 = -h(1)/(abs(h(1))^2+abs(h(2))^2);

y(2) = conj(y(2));

Y_eq = zeros(size(y));
Y_eq(1) = h11.*y(1) + h21.*y(2);
Y_eq(2) = h12.*y(1) + h22.*y(2);

Y_mod = qamdemod(Y_eq,4,'unitaveragepower',true,'outputtype','bit')';

received_sig = dec2bin(Y_mod);

errors = 0;
for k = [1:length(received_sig)]
    if received_sig(k) ~= data(k)
        errors = errors + 1;
    end
end
