clc;clear all;close all;
% generate data
N = 1024; %1000 bits
symbols = 2;
n_fft = 64;

data = dec2bin(round(rand(1,N)));
% modulate data using QPSK
data_sym = reshape(data,N/symbols,symbols);
data_dec = bin2dec(data_sym); 
mod_data = qammod(data_dec,4,'unitaveragepower',true);


% alamouti STBC
% [x1 -x2*; x2 x1*]
STBC = zeros(2,N/symbols);
STBC(:,1:2:end) = reshape(mod_data,2,N/(symbols*2));                        % creating [x1 ; x2]
STBC(:,2:2:end) = [-1;1].*flipud(reshape(conj(mod_data),2,N/(symbols*2)));  % creating [-x2* ; x1*]

% STBC_block1 = reshape(STBC(1,:),n_fft,length(STBC)/n_fft);
% STBC_block2 = reshape(STBC(2,:),n_fft,length(STBC)/n_fft);
% 
% % IFFT
% Tx1 = ifft(STBC_block1);
% Tx2 = ifft(STBC_block2);
% 
% % parallel to serial
% Tx1 = Tx1(:);
% Tx2 = Tx2(:);

% channel 
% x + j*y           [h1 ; h2]
h = [0.75+1*j; 0.9+0.8*j];


% freq domain, multiply channel with information

% y1 = conv(Tx1,h(1),'same');
% y2 = conv(Tx2,h(2),'same');


Tx1 = (STBC(1,:));
Tx2 = (STBC(2,:));
y1 = conv(Tx1,h(1),'same');
y2 = conv(Tx2,h(2),'same');

% superposition of two received signals at receiver
y = y1 + y2;

% serial to parallel
% y_sub = reshape(y,n_fft,length(STBC)/n_fft);
Y_blocks = y;
% DFT 
% Y_blocks  = fft(y_sub);
% lets say we estimated the channel and we know what h1 and h2 are

h11 = conj(h(1))/(abs(h(1))^2+abs(h(2))^2);
h21 = h(2)/(abs(h(1))^2+abs(h(2))^2);

h12 = conj(h(2))/(abs(h(1))^2+abs(h(2))^2);
h22 = -h(1)/(abs(h(1))^2+abs(h(2))^2);

Y_eq = zeros(size(Y_blocks));
Y_eq(:,1:2:end) = h11.*Y_blocks(:,1:2:end) + h21.*Y_blocks(:,2:2:end);
Y_eq(:,2:2:end) = h21.*Y_blocks(:,1:2:end) + h22.*Y_blocks(:,2:2:end);

Y_mod = qamdemod(Y_eq,4,'unitaveragepower',true);

received_sym = reshape(Y_mod,(length(data_sym)),1);
received_sig = reshape(dec2bin(received_sym),length(data),1);

errors = 0;
for k = [1:length(received_sig)]
    if received_sig(k) ~= data(k)
        errors = errors + 1;
    end
end

% % Receiver
% yMod = kron(reshape(y1,2,N),ones(1,2)); % [y1 y1 ... ; y2 y2 ...]
% yMod(2,:) = conj(yMod(2,:)); % [y1 y1 ... ; y2* y2*...]
% 
% % forming the equalization matrix
% hEq = zeros(2,N);
% hEq(:,[1:2:end]) = reshape(h,2,N/2); % [h1 0 ... ; h2 0...]
% hEq(:,[2:2:end]) = kron(ones(1,N/2),[1;-1]).*flipud(reshape(h,2,N/2)); % [h1 h2 ... ; h2 -h1 ...]
% hEq(1,:) = conj(hEq(1,:)); %  [h1* h2* ... ; h2 -h1 .... ]
% hEqPower = sum(hEq.*conj(hEq),1);
% 
% yHat = sum(hEq.*yMod,1)./hEqPower; % [h1*y1 + h2y2*, h2*y1 -h1y2*, ... ]
% yHat(2:2:end) = conj(yHat(2:2:end));
% 
% % receiver - hard decision decoding
% ipHat = real(yHat)>0;

% counting the errors
% nErr = size(find([data- ipHat]),2);