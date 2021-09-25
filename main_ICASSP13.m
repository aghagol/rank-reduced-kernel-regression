clc; clear;
% addpath(genpath('.\'));

srcdir = 'C:\Users\agahgol1\Desktop\SHARP\codes\samples_ICASSP13\';

for i = 1:6
    imname = ['patch', num2str(i), '.png'];
    y = imread([srcdir imname]);
    z = rrqr1(y);
    imwrite(z,[srcdir, imname(1:end-4) '_RRQR.png'])
end

