function [Ke,Fe,Re,Ve] = velocity(k)

if nargin==0
    Ke = [4 3 1 1 1 0];
    return
end

global U
global ELEMENT NODE DOF_ELEMENT
global BOUNDARY_DATA
global Flow

Ke = zeros(6);
Fe = zeros(6, 1);
Re = zeros(6);
Ve = zeros(6, 1);

% Integration points on [-1;+1]
N_GAUSS_POINT = 7;
a0 = 0;
a1 = 0.405845151377397;
a2 = 0.741531185599394;
a3 = 0.949107912342759;
c0 = 0.417959183673469;
c1 = 0.381830050505119;
c2 = 0.279705391489277;
c3 = 0.129484966168870;
GAUSS_POINT = [a0 -a1 a1 -a2 a2 -a3 a3];
GAUSS_WEIGHT = [c0 c1 c1 c2 c2 c3 c3];

% Node coordinates
xn = NODE(1:2, ELEMENT(1:6, k))';

% Current solution
phi = U(DOF_ELEMENT(1:6,k));

% Imposed normal velocity
V = BOUNDARY_DATA(1:3,k);

for n=1:N_GAUSS_POINT
    % Use the neighboring T6 element to calculate the velocity
    % Scale the coordinates because we use T6 shape functions defined on
    % the triangle between the vertices (0,0), (0,1) and (1,0).
    u = (GAUSS_POINT(n)+1)/2;
    v = 0;
    w = 1-u-v;
    J = [ 1-4*w 4*(w-u) 4*u-1 4*v 0 -4*v ; ...
        1-4*w -4*u 0 4*u 4*v-1 4*(w-v) ]*xn;
    GRAD = J\[ 1-4*w 4*(w-u) 4*u-1 4*v 0 -4*v ; ...
        1-4*w -4*u 0 4*u 4*v-1 4*(w-v) ];
    dphi = GRAD*phi;
    v2 = dphi(1)^2 + dphi(2)^2;
    rho = Flow.rho_inf*( 1 + (Flow.gamma-1)/2*(Flow.v_inf^2-v2)/Flow.c_inf^2 )^(1/(Flow.gamma-1));
    c2 = Flow.c_inf^2*( 1 + (Flow.gamma-1)/2*(Flow.v_inf^2-v2)/Flow.c_inf^2 );
    if v2>c2
        error('The flow is supersonic: the calculation cannot proceed (c=%f, v=%f).', sqrt(c2), sqrt(v2));
    end
    
    % L3 shape functions on [-1;+1]
    u = GAUSS_POINT(n);
    L = [-u*(1-u)/2 (1-u^2) u*(1+u)/2];
    Lu = [-1+2*u -4*u 1+2*u]/2;
    xs = Lu*xn(1:3, 1);
    ys = Lu*xn(1:3, 2);
    Js = sqrt(xs^2+ys^2);
    GRAD = Lu/Js;
    % Tangential vector
    tx = xs/Js;
    ty = ys/Js;
    % Gradients of the current solution in the tangential and normal
    % directions
    dphi_t = dphi(1)*tx + dphi(2)*ty;
    dphi_n = dphi(1)*ty - dphi(2)*tx;
    
    r = L*xn(1:3,2);
    
    Ke(1:3,1:3) = Ke(1:3,1:3) + GAUSS_WEIGHT(n)*( rho/c2*(L*V)*dphi_t*L'*GRAD )*r*Js;
    Fe(1:3) = Fe(1:3) + GAUSS_WEIGHT(n)*( -rho*L'*(L*V) + rho/c2*(L*V)*dphi_n*((L*V)-dphi_n)*L' )*r*Js;
end
