clc;clear all;close all;
% generate data
N = 4; 
symbols = 2;
%n_fft = 64;

data = round(rand(N,1));
% = [ 0 1 1 1]';
% dex = reshape(dec2bin(data'),2,2)';
% dec1 = bin2dec(dex);
% modulate data using QPSK

mod_data = qammod(data,4,'unitaveragepower',true,'inputtype','bit');          % 4QAM


% alamouti STBC
% [x1 -x2*; x2 x1*]                     
 

STBC = zeros(length(mod_data),1);
for i = 0:(length(mod_data)/2-1)
    STBC(4*i+1) = mod_data(2*i+1);              % creating [x1 ; x2]
    STBC(4*i+2) = mod_data(2*i+2);
    STBC(4*i+3) = -conj(mod_data(2*i+2));       % creating [-x2* ; x1*]
    STBC(4*i+4) = conj(mod_data(2*i+1));            
end
STBC = reshape(STBC,length(mod_data),2);

% channel 
% x + j*y           [h11 h12 ; h21 h22]
h = [randn()+randn()*1i, randn()+randn()*1i; randn()+randn()*1i, randn()+randn()*1i];
%h = [1 1; 1 1];

Tx1 = ifft(STBC(1,:));                   %[x1 -x2*]
Tx2 = ifft(STBC(2,:));                   %[x2 x1*]

y11 = Tx1*h(1,1);            % y time, receiver 
y21 = Tx2*h(2,1);

y12 = Tx1*h(1,2);
y22 = Tx2*h(2,2);

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
% assume perfect channel estimation
H = [ conj(h(1,1)) , conj(h(1,2)), h(2,1), h(2,2) ; ...             % (H^herm * H)^-1*H^herm 
    conj(h(2,1)), conj(h(2,2)), -h(1,1) , -h(1,2)]./ abs_h;        % pseudo inverse

Y1(2) = conj(Y1(2));
Y2(2) = conj(Y2(2));

%% this is the section that is different
% combined the expressions like the matrix
X_hat1 = H(1,1)*Y1(1) + H(1,2)*Y2(1) + H(1,3)*Y1(2) + H(1,4)*Y2(2);       
X_hat2 = H(2,1)*Y1(1) + H(2,2)*Y2(1) + H(2,3)*Y1(2) + H(2,4)*Y2(2);  

X(1) = X_hat1;
X(2) = X_hat2;

X_demod = qamdemod(X,4,'unitaveragepower',true,'outputtype','bit')';
output = reshape(X_demod.',size(data));


errors = 0;
for k = [1:N]
    if output(k) ~= data(k)
        errors = errors + 1;
    end
end
