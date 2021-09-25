function z = rrqr_vectorized(f)
g = double(rgb2gray(f));
f = double(f);
z = f;
q = zeros(size(f));
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
D1 = [xk(:) yk(:)] * [u1(:)'    ;u2(:)'];
D2 = [xk(:) yk(:)] * [-u2(:)'   ;u1(:)'];
W = exp(-.5 /sigma2 * bsxfun(@rdivide, bsxfun(@times,D1,e1(:)')+bsxfun(@times,D2,e2(:)'),alpha(:)'));
X11 = sum(W, 1);
X12 = sum(W .* D1, 1);
X22 = sum(W .* D1.^2, 1);
determinan = X11 .* X22 - X12.^2;
R11 =   X22 ./ determinan;
R12 =  -X12 ./ determinan;
% R22 =   X11 ./ determinan;
Xg1 = bsxfun(@times, W, R11) + bsxfun(@times, W .* D1, R12);
% Xg2 = R12 .* W + R22 .* W .* D1;
FR = im2col(padarray(f(:,:,1), [(pk-1)/2 (pk-1)/2], 'symmetric'), [pk pk], 'sliding');
FG = im2col(padarray(f(:,:,2), [(pk-1)/2 (pk-1)/2], 'symmetric'), [pk pk], 'sliding');
FB = im2col(padarray(f(:,:,3), [(pk-1)/2 (pk-1)/2], 'symmetric'), [pk pk], 'sliding');
ZR = sum(Xg1 .* FR, 1);
ZG = sum(Xg1 .* FG, 1);
ZB = sum(Xg1 .* FB, 1);
% QR = sum(Xg2 .* FR, 1);
% QG = sum(Xg2 .* FG, 1);
% QB = sum(Xg2 .* FB, 1);
z(:,:,1) = col2im(ZR, [pk pk], size(f(:,:,1)) + (pk-1), 'sliding');
z(:,:,2) = col2im(ZG, [pk pk], size(f(:,:,1)) + (pk-1), 'sliding');
z(:,:,3) = col2im(ZB, [pk pk], size(f(:,:,1)) + (pk-1), 'sliding');
% q(:,:,1) = col2im(QR, [pk pk], size(f(:,:,1)), 'sliding');
% q(:,:,2) = col2im(QG, [pk pk], size(f(:,:,1)), 'sliding');
% q(:,:,3) = col2im(QB, [pk pk], size(f(:,:,1)), 'sliding');

% for i = 1+(pk-1)/2 : size(g,1)-(pk-1)/2
% for j = 1+(pk-1)/2 : size(g,2)-(pk-1)/2
%     if alpha(i,j) > thresh
%         rb = f(i-(pk-1)/2:i+(pk-1)/2,j-(pk-1)/2:j+(pk-1)/2,1);
%         gb = f(i-(pk-1)/2:i+(pk-1)/2,j-(pk-1)/2:j+(pk-1)/2,2);
%         bb = f(i-(pk-1)/2:i+(pk-1)/2,j-(pk-1)/2:j+(pk-1)/2,3);
%         d1 = (xk * u1(i,j) + yk * u2(i,j)).^2;
%         d2 = (yk * u1(i,j) - xk * u2(i,j)).^2;
%         Wk = exp(-.5 / sigma2 / alpha(i,j) * (e1(i,j)*d1 + e2(i,j)*d2));
%         Ck = [ones(pk^2,1) d1(:)];
%         WCk = [Wk(:) d1(:).*Wk(:)];
%         Xg = (Ck' * WCk) \ WCk';
%         z(i,j,1) = Xg(1,:) * rb(:);
%         z(i,j,2) = Xg(1,:) * gb(:);
%         z(i,j,3) = Xg(1,:) * bb(:);
%         q(i,j,1) = Xg(2,:) * rb(:);
%         q(i,j,2) = Xg(2,:) * gb(:);
%         q(i,j,3) = Xg(2,:) * bb(:);
%     end
% end
% end
z = uint8(z - q.^5 .* exp(-.5 * q.^2 / 70) * .00008);
end