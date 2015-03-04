function [Ke,Fe,Re,Ve] = element(k)
% Compressible potential steady flow (3D)

if nargin==0
    Ke = [1 4 1 1 1 1 0];
    return
end

global Flow
global U
global ELEMENT NODE DOF_ELEMENT

% Integration points on the tetrahedron with vertices 
% (0,0,0), (0,1,0), (1,0,0) and (0,0,1)
N_GAUSS_POINT = 4;
a = (5-sqrt(5))/20; b = (5+3*sqrt(5))/20;
GAUSS_POINT = [ a a a b; a a b a; a b a a ];
GAUSS_WEIGHT = [1/24 1/24 1/24 1/24];

% Node coordinates
x_node = NODE(1,ELEMENT(1:4,k)).';
y_node = NODE(2,ELEMENT(1:4,k)).';
z_node = NODE(3,ELEMENT(1:4,k)).';

% Current solution
phi = U(DOF_ELEMENT(1:4,k));

Ke = zeros(4);
Fe = zeros(4, 1);
Re = zeros(4);
Ve = zeros(4, 1);

% Assembly
J = [-1 1 0 0 ; -1 0 1 0; -1 0 0 1]*[x_node y_node z_node];
detJ = det(J);
GRAD = J\[-1 1 0 0 ; -1 0 1 0; -1 0 0 1];

for n=1:N_GAUSS_POINT
    %u = GAUSS_POINT(1,n);
    %v = GAUSS_POINT(2,n);
    %w = GAUSS_POINT(3,n);
    %************************
    % N = [(1-u-v-w)  u  v  w];
    
    dphi = GRAD*phi;
    v2 = dphi(1)^2 + dphi(2)^2 + dphi(3)^2;
    rho = Flow.rho_inf*( 1 + (Flow.gamma-1)/2*(Flow.v_inf^2-v2)/Flow.c_inf^2 )^(1/(Flow.gamma-1));
    c2 = Flow.c_inf^2*( 1 + (Flow.gamma-1)/2*(Flow.v_inf^2-v2)/Flow.c_inf^2 );
    if v2>c2
        error('The flow is supersonic: the calculation cannot proceed (c=%f, v=%f).', sqrt(c2), sqrt(v2));
    end
    Ke = Ke + GAUSS_WEIGHT(n)*( rho/c2*(dphi'*GRAD)'*(dphi'*GRAD) - rho*(GRAD'*GRAD) )*detJ;
    Fe = Fe + GAUSS_WEIGHT(n)*( rho*GRAD'*dphi )*detJ;
end


