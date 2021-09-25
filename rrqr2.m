function [z q] = rrqr2(pk,g,alpha,thresh,f,xk,yk,u1,u2,sigma2,e1,e2)
    z = f;
    q = zeros(size(f));
    for i = 1+(pk-1)/2 : size(g,1)-(pk-1)/2
    for j = 1+(pk-1)/2 : size(g,2)-(pk-1)/2
        if alpha(i,j) > thresh
            rb = f(i-(pk-1)/2:i+(pk-1)/2,j-(pk-1)/2:j+(pk-1)/2,1);
            gb = f(i-(pk-1)/2:i+(pk-1)/2,j-(pk-1)/2:j+(pk-1)/2,2);
            bb = f(i-(pk-1)/2:i+(pk-1)/2,j-(pk-1)/2:j+(pk-1)/2,3);
            d1 = (xk * u1(i,j) + yk * u2(i,j)).^2;
            d2 = (yk * u1(i,j) - xk * u2(i,j)).^2;
            Wk = exp(-.5 / sigma2 / alpha(i,j) * (e1(i,j)*d1 + e2(i,j)*d2));
            Ck = [ones(pk^2,1) d1(:)];
            WCk = [Wk(:) d1(:).*Wk(:)];
            Xg = (Ck' * WCk) \ WCk';
            z(i,j,1) = Xg(1,:) * rb(:);
            z(i,j,2) = Xg(1,:) * gb(:);
            z(i,j,3) = Xg(1,:) * bb(:);
            q(i,j,1) = Xg(2,:) * rb(:);
            q(i,j,2) = Xg(2,:) * gb(:);
            q(i,j,3) = Xg(2,:) * bb(:);
        end
    end
    end
end