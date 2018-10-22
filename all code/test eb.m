clc;clear all;close all;

eb = [0:1:30];

ber = berawgn(eb, '4qam');

plot(eb,ber)