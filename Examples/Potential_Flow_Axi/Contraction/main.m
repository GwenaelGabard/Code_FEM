% Solve the steady, potential, compressible flow through an axi-symmetric duct

% Add the Core directory to the path
addpath(fullfile('..','..','..','..','Code_FEM','Core'));

global LIBRARY_PATH %#ok<NUSED>
MANAGE_PATH

% Properties of the flow at infinity
global Flow
Flow.gamma = 1.4;
Flow.rho_inf = 1.2;
Flow.c_inf = 340;
Flow.outlet_velocity = Flow.c_inf*0.25;
Flow.v_inf = Flow.outlet_velocity;

% Create a T3 mesh
% [node, edge, element] = initmesh('geometry', 'Jiggle', 'mean', 'JiggleIter', 20, 'Hmax', 0.08);
% or load an example of mesh
load mesh_example

% Convert the mesh to T6 elements
[node, edge, element] = convert_T3_T6(node, edge, element, @geometry);

% Plot the mesh
figure;
PLOT_FEM(node(1:2, :), element(1:6, :));
axis equal; axis tight; box on;

INIT_GEOMETRY

ADD_DOMAIN(node, edge, element);
 
BUILD_MODEL
BUILD_DOF

% The boundary elements need to have access to the T6 elements
ADD_BORDER_DOFS(1, 1:8);


%% Calculate the incompressible flow

% The current solution will be used by the element functions
global U

% Initial solution
U = zeros(N_DOF,1);

% Instead of using a separate model for the incompressible flow the trick
% is to set the reference sound speed to infinity. All the terms associated
% with compressibility will naturally drop from the formulation.
c_inf_backup = Flow.c_inf;
Flow.c_inf = inf;

% Build and solve the linear system
BUILD_SYSTEM;
SOLVE_SYSTEM;

% Restore the actual value of the reference sound speed
Flow.c_inf = c_inf_backup;
clear c_inf_backup


%% Solving non-linear problem for the compressible flow

% Maximum number of iterations for each step
N_MAX_ITER = 100;
% Convergence criteria for the iterations
ERROR = 1.0e-12;

% We use a basic Newton-Raphson method to solve the non-linear equations
for iter=1:N_MAX_ITER
    fprintf('\nITERATION %i\n', iter);
    % Build and solve the linear system for this iteration
    BUILD_SYSTEM;
    Uold = U;
    SOLVE_SYSTEM;
    % The norm of the delta
    e = norm(U);
    % Update the solution
    U = Uold+U;
    % Relative change in solution
    e = e/norm(U);
    % Test for convergence
    fprintf('ITERATION ERROR: %e\n',e);
    if e < ERROR
        break
    end
end

% Plot the velocity potential
figure
PLOT_FEM(node(1:2,:), element(1:6,:), full(U))
axis equal; axis tight; box on;
colorbar;
title('velocity potential');


%% Calculate the velocity and other flow properties
get_flow;

% Plot the properties of the flow
figure
PLOT_FEM(node(1:2,:), element(1:6,:), Vx)
axis equal; axis tight; box on;
colorbar;
title('vx');

figure
PLOT_FEM(node(1:2,:), element(1:6,:), Vy)
axis equal; axis tight; box on;
colorbar;
title('vy');

figure
PLOT_FEM(node(1:2,:), element(1:6,:), V./c)
axis equal; axis tight; box on;
colorbar;
title('Mach number');

figure
PLOT_FEM(node(1:2,:), element(1:6,:), rho)
axis equal; axis tight; box on;
colorbar;
title('density');

figure
PLOT_FEM(node(1:2,:), element(1:6,:), c)
axis equal; axis tight; box on;
colorbar;
title('sound speed');
