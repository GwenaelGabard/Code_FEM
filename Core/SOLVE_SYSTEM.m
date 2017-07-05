% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% Copyright (c) 2014 Code_FEM developers


% Solve the linear system and implement the linear constraints

fprintf('SOLVE_SYSTEM: ');
% Find the DOFs that are constrained
L = find(diag(R));
fprintf('%i constraints', length(L));
% Find the DOFs that are not constrained
Lbar = find(diag(R)==0);

% Change of variables to impose the linear constraints
R = R + sparse(Lbar, Lbar, ones(1, length(Lbar)), N_DOF, N_DOF);
R = inv(R);
K = R.'*K*R;
F = R.'*F;

% Modify the linear system 
% K(L,:) = 0;% K((L-1)*N_DOF+L) = 1;% F(L) = V(L); (original method - unsymmetric + slow)
F = F-(K(:,L)*V(L)); K(L,:) = 0; K(:,L) = 0; 
diagK = diag(K); diagK(L)=1; K = spdiags(diagK,0,K); 
F(L,1:N_RHS) = V(L,1:N_RHS);

% Direct resolution of the linear system
global U
U = K\F;

% Convert back to the original variables
U = R*U;
fprintf(', %i DOFs\n',N_DOF);

% Clean up the workspace
clear L Lbar k
