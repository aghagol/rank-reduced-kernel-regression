clc; clear;
addpath(genpath('.\'));

srcdir = 'C:\Users\agahgol1\Desktop\SHARP\codes\samples_ICASSP13\';
% srcfiles = dir([srcdir, '*.png']);

ksize1  = 5;
ksize2  = 9;
wsize   = 9;
lambda  = 1;
alpha   = .4;
r       = 1;
h1      = .4;
h2      = 3;

for i = 1:6 %numel(srcfiles)
%     imname = srcfiles(i).name;
    imname = ['patch', num2str(i), '.png'];
    y = imread([srcdir imname]);
    g = double(rgb2gray(y));
    y = double(y);
    z = y*0;
    [gg, gx, gy] = ckr2_regular(g, h1, r, ksize1);
    C = steering(gx, gy, ones(size(g)), wsize, lambda, alpha);
    for k = 1:3
        zz = skr0_regular(y(:,:,k), h2, C, r, ksize2);
        z(:,:,k) = zz;
    end
    imwrite(uint8(z),[srcdir, imname(1:end-4) '_SKR0.png'])
end

