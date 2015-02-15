% Calculate the flow properties at the nodes of the T6 mesh

% Coordinates of the points
UInt = [0 0.5 1 0.5 0 0];
VInt = [0 0 0 0.5 1 0.5];

% Shape functions and their derivatives
WInt = 1 - UInt - VInt;
N = [ (2*WInt-1).*WInt ; 4*UInt.*WInt ; (2*UInt-1).*UInt ; 4*UInt.*VInt ; (2*VInt-1).*VInt ; 4*VInt.*WInt ];
dNdu = [ 1-4*WInt ; 4*(WInt-UInt) ; 4*UInt-1 ; 4*VInt ; 0*UInt ; -4*VInt ];
dNdv = [ 1-4*WInt ; -4*UInt ; 0*UInt ; 4*UInt ; 4*VInt-1 ; 4*(WInt-VInt) ];

% In each element we calculate the gradient of the velocity potential.
% We end up with several values of the velocity since the nodes are shared
% by several elements. The average of these values is calculated.
ElemList = find(ELEMENT_DOMAIN(2,:)==0);
Vx = zeros(1, size(NODE, 2));
Vy = Vx;
NV = Vx;
for k=ElemList
    nn = ELEMENT(1:6,k);
    xn = NODE(1:2, nn).';
    Ue = full(U(DOF_ELEMENT(1:N_DOF_ELEMENT(k),k)).').';
    for n=1:length(nn)
        GRAD = [dNdu(:,n) dNdv(:,n)]';
        J = GRAD*xn;
        GRAD = J\GRAD;
        dphi = GRAD*Ue;
        Vx(nn(n)) = Vx(nn(n)) + dphi(1);
        Vy(nn(n)) = Vy(nn(n)) + dphi(2);
    end
    NV(nn) = NV(nn) +1;
end
% Average the values
Vx = Vx./NV;
Vy = Vy./NV;

% Velocity normal
V2 = Vx.^2 + Vy.^2;
V = sqrt(V2);

% Density
rho = Flow.rho_inf*( 1 + (Flow.gamma-1)/2*(Flow.v_inf^2-V2)/Flow.c_inf^2 ).^(1/(Flow.gamma-1));

% Sound speed
c2 = Flow.c_inf^2*( 1 + (Flow.gamma-1)/2*(Flow.v_inf^2-V2)/Flow.c_inf^2 );
c = sqrt(c2);
