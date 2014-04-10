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

function [Ke,Fe,Re,Ve] = element(k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Convected Helmholtz Equation
% Notation: exp(-i \omega t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear triangular element
% 2 nodes
% 3 DOFs: p1, p2 and p3
% Data in ELEMENT_DATA(:,k):
% None
% Data in NODE_DATA(:,k):
% rho0, c0, u0, v0, s

if nargin==0
    % Return the properties of the element
    % Ke should contain:
    % the priority of the element
    % the number of nodes
    % the number of DOFs for node 1
    % the number of DOFs for node 2
    % ...
    % the number of internal DOFs of the element
   Ke = [1 4 1 1 1 1 0];
   return
end

global omega
global ELEMENT NODE NODE_DATA

% Sound speed at the nodes
C0 = NODE_DATA(2,ELEMENT(1,k));
% Mean flow velocity at the nodes
U0 = NODE_DATA(3:4,ELEMENT(1:3,k));
% Source term at the nodes
S = NODE_DATA(5,ELEMENT(1:3,k));

% Gauss points on the reference element
N_GAUSS_POINT = 4;
GAUSS_POINT = [1/sqrt(3) 1/sqrt(3) -1/sqrt(3) -1/sqrt(3); 1/sqrt(3) -1/sqrt(3) 1/sqrt(3) -1/sqrt(3) ];
GAUSS_WEIGHT = [1 1 1 1];

% Coordinates of the nodes
x_node = NODE(1,ELEMENT(1:4,k)).';
y_node = NODE(2,ELEMENT(1:4,k)).';

% Initialise the output
Ke = zeros(4,4);
Fe = zeros(4,1);
Re = zeros(4,4);
Ve = zeros(4,1);

% Integrate over the reference elements using the Gauss points
for n=1:N_GAUSS_POINT
   % Coordinates of the Gauss point on the reference elements
   u = GAUSS_POINT(1,n);
   v = GAUSS_POINT(2,n);
   % Jacobian to the physical element
   J = [v-1 1-v 1+v -1-v ;u-1 -u-1 1+u 1-u]/4*[x_node y_node];
   detJ = det(J);

   % Shape functions for the solution
   N = [(1-u)*(1-v) (1+u)*(1-v) (1+u)*(1+v) (1-u)*(1+v)]/4;
   % Shape functions for the gradient of the solution
   GRAD = J\[v-1 1-v 1+v -1-v ;u-1 -u-1 1+u 1-u]/4;
   % Material derivative
    D0pDt = -1i*omega*N + (N*U0')*GRAD;
    
    % Add the contribution to the element matrix
    Ke = Ke + GAUSS_WEIGHT(n)*( C0^2*(GRAD'*GRAD) - D0pDt'*D0pDt )*detJ;
    % Add the contribution to the right-hand side
    Fe = Fe + GAUSS_WEIGHT(n)*( N.'*N*S.' )*detJ;
end

