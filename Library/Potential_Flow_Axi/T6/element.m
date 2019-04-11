function [Ke,Fe,Re,Ve] = element(k)
% Compressible potential steady flow
% Axisymmetric flow
% Cylindrical coordinates: z, r


if nargin==0
    Ke = [1 6 1 1 1 1 1 1 0];
    return
end

global Flow
global U
global ELEMENT NODE DOF_ELEMENT

% Integration points on the triangle with vertices (0,0), (0,1) and (1,0)
N_GAUSS_POINT = 7;
GAUSS_POINT = [ 0.797426985353087   0.101286507323456 ; ...
    0.101286507323456   0.797426985353087 ; ...
    0.101286507323456   0.101286507323456 ; ...
    0.059715871789770   0.470142064105115 ; ...
    0.470142064105115   0.059715871789770 ; ...
    0.470142064105115   0.470142064105115 ; ...
    0.333333333333333   0.333333333333333 ]';
GAUSS_WEIGHT = [   0.125939180544827
    0.125939180544827 ; ...
    0.125939180544827 ; ...
    0.132394152788506 ; ...
    0.132394152788506 ; ...
    0.132394152788506 ; ...
    0.225000000000000]'/2;

% Node coordinates
xn = NODE(1:2,ELEMENT(1:6,k))';

% Current solution
phi = U(DOF_ELEMENT(1:6,k));

Ke = zeros(6);
Fe = zeros(6, 1);

for n=1:N_GAUSS_POINT
    u = GAUSS_POINT(1,n);
    v = GAUSS_POINT(2,n);
    w = 1-u-v;
    J = [ 1-4*w 4*(w-u) 4*u-1 4*v 0 -4*v ; ...
        1-4*w -4*u 0 4*u 4*v-1 4*(w-v) ]*xn;
    detJ = det(J);
    GRAD = J\[ 1-4*w 4*(w-u) 4*u-1 4*v 0 -4*v ; ...
        1-4*w -4*u 0 4*u 4*v-1 4*(w-v) ];
    N = [ (2*w-1)*w 4*u*w (2*u-1)*u 4*u*v (2*v-1)*v 4*v*w ];
    r = N*xn(:,2);
    
    dphi = GRAD*phi;
    v2 = dphi(1)^2 + dphi(2)^2;
    rho = Flow.rho_inf*( 1 + (Flow.gamma-1)/2*(Flow.v_inf^2-v2)/Flow.c_inf^2 )^(1/(Flow.gamma-1));
    c2 = Flow.c_inf^2*( 1 + (Flow.gamma-1)/2*(Flow.v_inf^2-v2)/Flow.c_inf^2 );
    if v2>c2
        warning('The flow is supersonic: the calculation cannot proceed (c=%f, v=%f).', sqrt(c2), sqrt(v2));
    end
    Ke = Ke + GAUSS_WEIGHT(n)*( rho/c2*(dphi'*GRAD)'*(dphi'*GRAD) - rho*(GRAD'*GRAD) )*r*detJ;
    Fe = Fe + GAUSS_WEIGHT(n)*( rho*GRAD'*dphi )*r*detJ;
end

Re = zeros(6);
Ve = zeros(6, 1);
