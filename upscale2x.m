clc;
srcdir = 'C:\Users\agahgol1\Desktop\SHARP\Groundtruth\';
dstdir = 'C:\Users\agahgol1\Desktop\SHARP\Upscaled2X\';
seqtag = 'flower';
srcdir = [srcdir seqtag '\'];
dstdir = [dstdir seqtag '\'];
if ~isdir(dstdir)
    mkdir(dstdir);
end
imlist = dir([srcdir '*.png']);
for i = 1:numel(imlist)
    x = double(imread([srcdir imlist(i).name]));
    y = dnscale(x);
    z = upscale(y);
    x(1:end-1,1:end-1,:) = z;
    imwrite(uint8(x), [dstdir 'BL_' imlist(i).name]);
end