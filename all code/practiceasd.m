clc;clear all;close all;

data = ones(100,1);

window = 10;
delay = 2;
for i = 0:length(data)/window-2
    data(1+window+window*i:1+window+delay+window*i) = data(1+window+window*i:1+window+delay+window*i) ...
        + data(window-delay+window*i:window+window*i)
end

data(1:window) = [zeros(delay+1,1);data(1:window-delay-1)];       