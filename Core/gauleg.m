function [x,w] = gauleg(x1,x2,n)
% Gauss legendre quadrature of order n between x1 and x2
m = (n+1)/2;
xm = 0.5*(x2+x1);
xl = 0.5*(x2-x1);

EPS = 3.0E-14; 

x = zeros(1,n);
w = zeros(1,n);
for i=1:m
    z = cos(pi*(i-0.25)/(n+0.5));
    z1 = z +10;
    while abs(z-z1)>EPS
        p1 = 1;
        p2 = 0;
        for j =1:n
            p3 = p2;
            p2 = p1;
            p1 = ((2*j-1)*z*p2-(j-1)*p3)/j;
        end
        pp = n*(z*p1-p2)/(z*z-1);
        z1 = z;
        z = z1-p1/pp;
    end
    x(i) = xm-xl*z;
    x(n+1-i) = xm+xl*z;
    w(i) = 2*xl/((1-z*z)*pp*pp);
    w(n+1-i) = w(i);
end

return
