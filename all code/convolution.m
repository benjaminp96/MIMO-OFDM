clc;close all

vector1 = [1 2 3 4 5];
vector2 = [6 7 8 9 10];

channel1 = [ .5+.75*1j];
channel2 = [ .75+.5*1j];
conv1 = conv(vector1,channel1);
conv2 = conv(vector2,channel2);
vector = [vector1;vector2];
channel = [ channel1;channel2];
conv3 = conv(vector1,channel);