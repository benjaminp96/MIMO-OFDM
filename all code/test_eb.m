clc;clear all;close all;

eb = [0:1:20];
div = [1:1:20];
ber1 = berawgn(eb, 'qam',4);
ber2 = berawgn(eb, 'qam',16);
ber3 = berawgn(eb, 'qam',64);

semilogy(eb,ber1)
hold on
semilogy(eb,ber2)
semilogy(eb,ber3)

ber1 = berfading(eb, 'qam',4,div);
ber2 = berfading(eb, 'qam',16,div);
ber3 = berfading(eb, 'qam',64,div);
figure


semilogy(eb,ber1)
hold on
semilogy(eb,ber2)
semilogy(eb,ber3)