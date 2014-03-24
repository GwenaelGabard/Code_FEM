% To find the solver
addpath('../../../Code_FEM/Core');
% to do add an environment variable with the path_name
MANAGE_PATH

global omega

omega = 10; 
C0 = 1;
RHO0 = 1;
Length = 1;
n_elem = 10000;

INIT_GEOMETRY
[N_DIM,N_DOMAIN,N_NODE,NODE,NODE_ID,NODE_DOMAIN,...
N_ELEMENT,ELEMENT,ELEMENT_ID,ELEMENT_DOMAIN,N_NODE_ELEMENT] = MESH1D(Length,n_elem);
BUILD_MODEL
BUILD_DOF
BUILD_SYSTEM

% admitance BC
K(end,end) = K(end,end) + 1i*omega;
SOLVE_SYSTEM
U = full(U);

% Analytical solution
Uexact = exp(-1i*omega*real(NODE)).';
%
figure;
plot(real(NODE),real(U),'k.'); hold on;
plot(real(NODE),real(Uexact),'k-','Markersize',5); hold on;
set(gca,'Fontsize',14);
xlabel('x');
ylabel('real part of pressure');

NodalError = 100*norm(U(1:n_elem+1)-Uexact(1:n_elem+1))./norm(Uexact(1:n_elem+1));
fprintf('Erreur = %g pc\n',NodalError);
