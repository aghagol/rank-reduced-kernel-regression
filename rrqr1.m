function z = rrqr1(f)
g = double(rgb2gray(f));
f = double(f);
%% set parameters
pg = 5;
ps = 9;
pk = 9;
thresh = 15;
sigma2 = 3^2;
%% compute gradient field for luminance
[xg yg] = meshgrid(-(pg-1)/2:(pg-1)/2);
Wg = exp(-.5 / .4^2 * (xg.^2 + yg.^2));
Cg = [xg(:) yg(:)];
WCg = [xg(:).*Wg(:) yg(:).*Wg(:)];
Xg = (Cg' * WCg) \ WCg';
gx = imfilter(g, reshape(Xg(1,:),[pg pg]), 'symmetric', 'same', 'corr');
gy = imfilter(g, reshape(Xg(2,:),[pg pg]), 'symmetric', 'same', 'corr');
%% compute structure tensor field
[xs ys] = meshgrid(-(ps-1)/2:(ps-1)/2);
Ws = exp(-.5 / 3^2 * (xs.^2 + ys.^2));
Ws = Ws ./ sum(Ws(:));
s1 = imfilter(gx .* gx, Ws, 'symmetric', 'same', 'corr');
s2 = imfilter(gx .* gy, Ws, 'symmetric', 'same', 'corr');
s4 = imfilter(gy .* gy, Ws, 'symmetric', 'same', 'corr');
%% compute structure spectrum for all pixels
delta_sq = sqrt((s1 - s4).^2 + 4 * s2.^2);
e1 = .5 * abs(s1+s4 + delta_sq);
e2 = .5 * abs(s1+s4 - delta_sq);
u2 = (e1 - s1) ./ s2;
temp = sqrt(1 + u2.^2);
u1 = 1  ./ temp;
u2 = u2 ./ temp;
%% compute kernel variance
alpha = ((e1 - e2)./(e1 + e2)).^2 .* (e1 .* e2).^.4;
e1 = min(e1, 10 * alpha * sigma2);
%% rank-reduced regression
[xk yk] = meshgrid(-(pk-1)/2:(pk-1)/2);
[z q] = mexrrqr2(pk,g,alpha,thresh,f,xk,yk,u1,u2,sigma2,e1,e2);
z = uint8(z - q.^5 .* exp(-.5 * q.^2 / 70) * .00008);
end