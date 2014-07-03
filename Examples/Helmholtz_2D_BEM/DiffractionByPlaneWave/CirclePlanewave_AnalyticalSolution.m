function [U] = CirclePlanewave_AnalyticalSolution(node,a)

global omega

% Analytical solution of the diffraction of a plane wave by a circle (zero
% incidence)

k = omega/340;
[T,R] = cart2pol(node(1,:),node(2,:));
Phiinc = exp(1i*k*node(1,:)); 
r = a;
Phida = - ((besselj(1,k*a))/besselh(1,1,k*a))*besselh(0,1,k*r)+0*length(T);

for temp = 1:10000
    dJa = 0.5*k*(besselj(temp-1,k*a)-besselj(temp+1,k*a));
    dYa = 0.5*k*(bessely(temp-1,k*a)-bessely(temp+1,k*a));
    dHa = dJa+1i*dYa;
    temp2 = - 2*1i^temp*cos(temp*T)*(dJa/dHa)*besselh(temp,1,k*r);
    if isnan(temp2)==0
        Phida = Phida + temp2;
    else
        break
    end
end

% Total field
U = (Phida + Phiinc).';
