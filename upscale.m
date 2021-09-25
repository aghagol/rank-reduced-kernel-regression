function z = upscale(y)
    z = zeros([2*size(y,1)-1 2*size(y,2)-1 size(y,3)]);
    for i = 1:size(y,3)
        z(:,:,i) = interp2(y(:,:,i));
    end
end