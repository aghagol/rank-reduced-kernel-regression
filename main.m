% clc; 
clear;
% dbstop if error
srcdir = 'C:\Users\agahgol1\Desktop\SHARP\Upscaled2X\';
dstdir = 'C:\Users\agahgol1\Desktop\SHARP\test\';
seqtag = 'NBA';
srcdir = [srcdir seqtag '\'];
dstdir = [dstdir seqtag '\'];
if ~isdir(dstdir)
    mkdir(dstdir);
end
imlist = dir([srcdir '*.png']);
for i = 1%:numel(imlist)
    iname = imlist(i).name;
    x = imread([srcdir iname]);
    tic
    y = rrqr1(x);
    toc
%     imwrite(y, [dstdir 'QR_' iname(4:end)]);
end
