function [Ke,Fe,Re,Ve] = velocity(k)

if nargin==0
    Ke = [4 3 1 1 1 0];
    return
end

global U
global ELEMENT NODE DOF_ELEMENT
global BOUNDARY_DATA
global Flow

Ke = zeros(4);
Fe = zeros(4,1);
Re = zeros(4);
Ve = zeros(4,1);

% Integration points on triangle with vertex (0,0), (0,1) and (1,0)
N_GAUSS_POINT = 3;
GAUSS_POINT = [1/2 0 1/2 ; 1/2 1/2 0];
GAUSS_WEIGHT = [1/6 1/6 1/6];
%
% Node coordinates (tetra4)
x_node = NODE(1,ELEMENT(1:4,k)).';
y_node = NODE(2,ELEMENT(1:4,k)).';
z_node = NODE(3,ELEMENT(1:4,k)).';

% Node coordinates (tria3)
x_point = NODE(1,ELEMENT(1:3,k)).';
y_point = NODE(2,ELEMENT(1:3,k)).';
z_point = NODE(3,ELEMENT(1:3,k)).';

% Geometry for integration (tria3) 
tu = [[-1 1 0]*x_point; [-1 1 0]*y_point; [-1 1 0]*z_point]; 
tv = [[-1 0 1]*x_point; [-1 0 1]*y_point; [-1 0 1]*z_point];
tw = [tu(2)*tv(3)-tu(3)*tv(2) ; tu(3)*tv(1)-tu(1)*tv(3) ; tu(1)*tv(2)-tu(2)*tv(1)];
Js = norm(tw); 
no = tw/norm(tw);

% Current solution
phi = U(DOF_ELEMENT(1:4,k));

% Imposed normal velocity
Vnodes = BOUNDARY_DATA(1:3,k);

for n=1:N_GAUSS_POINT
    % Use the neighboring T4 element to calculate the velocity
    u = GAUSS_POINT(1,n);
    v = GAUSS_POINT(2,n);
    w = 0;
    %
    N = [(1-u-v-w)  u  v  w];
    X_tetra = N*NODE(1:3,ELEMENT(1:4,k)).';
    % 
    J = [-1 1 0 0 ; -1 0 1 0; -1 0 0 1]*[x_node y_node z_node];
    GRAD = J\[-1 1 0 0 ; -1 0 1 0; -1 0 0 1];
    %
    dphi = (GRAD*phi).';
    v2 = dphi(1)^2 + dphi(2)^2 + dphi(3)^2 ;
    rho = Flow.rho_inf*( 1 + (Flow.gamma-1)/2*(Flow.v_inf^2-v2)/Flow.c_inf^2 )^(1/(Flow.gamma-1));
    c2 = Flow.c_inf^2*( 1 + (Flow.gamma-1)/2*(Flow.v_inf^2-v2)/Flow.c_inf^2 );
    if v2>c2
        error('The flow is supersonic: the calculation cannot proceed (c=%f, v=%f).', sqrt(c2), sqrt(v2));
    end
    
    % T3 shape functions
    u = GAUSS_POINT(1,n);
    v = GAUSS_POINT(2,n);
    %
    L = [1-u-v u v]; 
    X_tria = L*NODE(1:3,ELEMENT(1:3,k)).';
    if norm(X_tetra-X_tria)>100*eps
        error('Problem in velocity!')
    end

    % Gradients of the current solution in the normal direction
    dphi_n = dphi(1)*no(1) + dphi(2)*no(2) + dphi(3)*no(3);
    
    % normal gradient
    Gn = GRAD(1,:)*no(1) + GRAD(2,:)*no(2) + GRAD(3,:)*no(3);
    
    % local prescribed velocity 
    V = L*Vnodes;
    
    Ke(1:4,1:4) = Ke(1:4,1:4) + GAUSS_WEIGHT(n)*( rho/c2*V*N.'*(   dphi*GRAD - dphi_n*Gn  )    )*Js;
    Fe(1:4) = Fe(1:4) + GAUSS_WEIGHT(n)*( -rho*N'*V + rho/c2*N.'*V*dphi_n*(V-dphi_n) )*Js;
end