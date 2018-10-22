clc;clear all;close all;
% generate data
N = 4; %1000 bits
symbols = 2;
%n_fft = 64;

data = dec2bin(round(rand(N,1)));
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
% x + j*y           [h11 h12 ; h21 h22]
h = [randn()+randn()*1i, randn()+rand()*1i; randn()+rand()*1i, randn()+rand()*1i];
h = h./norm(h);


Tx1 = ifft(STBC(1,:));                   %[x1 -x2*]
Tx2 = ifft(STBC(2,:));                   %[x2 x1*]

y11 = conv(Tx1,h(1,1),'same');            % y time, receiver 
y21 = conv(Tx2,h(2,1),'same');

y12 = conv(Tx1,h(1,2),'same');
y22 = conv(Tx2,h(2,2),'same');

% y1 = Tx1*h(1);
% y2 = Tx2*h(2);

% superposition of two received signals at receiver
y1 = y11 + y21;                    % y(1) = h1*x1 + h2*x2 ; y(2) = -h2* x1* + h1 * x2*
y2 = y12 + y22;
% serial to parallel
% y_sub = reshape(y,n_fft,length(STBC)/n_fft);

% DFT 
% Y_blocks  = fft(y_sub);
% lets say we estimated the channel and we know what h1 and h2 are

Y1 = fft(y1);
Y2 = fft(y2);

abs_h = sum(sum(abs(h).^2));

H = ([ conj(h(1,1)) , conj(h(1,2)), h(2,1), h(2,2) ; ...
    conj(h(2,1)), conj(h(2,2)), -h(1,1) , -h(1,2)]./ abs_h);

Y1(2) = conj(Y1(2));
Y2(2) = conj(Y2(2));

Y_eq1 = zeros(size(Y1));
Y_eq1(1) = H(1,1).*Y1(1) + H(1,3).*Y1(2);       
Y_eq1(2) = H(2,1).*Y1(1) + H(2,3).*Y1(2);


Y_eq2 = zeros(size(Y2));
Y_eq2(1) = H(1,2).*Y2(1) + H(1,4).*Y2(2);
Y_eq2(2) = H(2,2).*Y2(1) + H(2,4).*Y2(2);

Y_mod1 = qamdemod(Y_eq1,4,'unitaveragepower',true,'outputtype','bit')';

received_sig1 = dec2bin(Y_mod1);

errors1 = 0;
for k = [1:N]
    if received_sig1(k) ~= data(k)
        errors1 = errors1 + 1;
    end
end
Y_mod2 = qamdemod(Y_eq2,4,'unitaveragepower',true,'outputtype','bit')';

received_sig2 = dec2bin(Y_mod2);

errors2 = 0;
for k = [1:N]
    if received_sig2(k) ~= data(k)
        errors2 = errors2 + 1;
    end
end