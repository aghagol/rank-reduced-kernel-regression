clc; clear;

z = [];
for i = 1:6
    imname = ['patch', num2str(i), '.png'];
    y = imread(imname);
    temp = y(6:end-5,6:end-5,:);
    
    imname = ['patch', num2str(i), '_SKR2.png'];
    y = imread(imname);
    temp = cat(1,temp,y(6:end-5,6:end-5,:));
    
    imname = ['patch', num2str(i), '_SKR0.png'];
    y = imread(imname);
    temp = cat(1,temp,y(6:end-5,6:end-5,:));
    
    imname = ['patch', num2str(i), '_RRQR.png'];
    y = imread(imname);
    temp = cat(1,temp,y(6:end-5,6:end-5,:));
    
    z = cat(2,z,temp);
end

imwrite(z,'merged.png')

imshow(z)
saveas(gca, 'merged.eps', 'psc2')