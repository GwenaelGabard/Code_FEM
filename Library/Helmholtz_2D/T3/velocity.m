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


function [Ke,Fe,Re,Ve] = velocity(k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Convected Helmholtz Equation
% Notation: exp(+i \omega t)
% Boundary condition: velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear boundary element associated to the T3 element
% 2 nodes
% 2 DOFs: p1, p2
% Data in ELEMENT(:,k) :
% 1,2: node numbers
% 3,4: coordinates along the boundary
% 5: number of this boundary
% 6: number of the neighboring element
% 7: number of the side in the neighboring element
% Data in ELEMENT_DATA(:,k):
% None
% Data in NODE_DATA(:,k):
% None
% Data in BOUNDARY_DATA(:,k):
% 1: impedance at node 1
% 2: impedance at node 2

if nargin==0
    % Return the properties of the element
    % Ke should contain:
    % the priority of the element
    % the number of nodes
    % the number of DOFs for node 1
    % the number of DOFs for node 2
    % ...
    % the number of internal DOFs of the element
    Ke = [5 2 1 1 0];
    return
end

global omega
global ELEMENT NODE NODE_DATA BOUNDARY_DATA

RHO = NODE_DATA(1,ELEMENT(1,k));
C0 = NODE_DATA(2,ELEMENT(1,k));

% Boundary condition values (velocity at nodes)
Vx1 = BOUNDARY_DATA(1,k);
Vy1 = BOUNDARY_DATA(2,k);
Vx2 = BOUNDARY_DATA(3,k);
Vy2 = BOUNDARY_DATA(4,k);

Ke = zeros(2,2);
Fe = zeros(2,1);
Re = zeros(2,2);
Ve = zeros(2,1);

N_GAUSS_POINT = 2;
GAUSS_POINT = [1 -1]/sqrt(3);
GAUSS_WEIGHT = [1 1];

x_point = NODE(1,ELEMENT(1:2,k))';
y_point = NODE(2,ELEMENT(1:2,k))';

% Compute the normal
nox = [-1 1]*y_point;
noy = [1 -1]*x_point;
no = [nox noy]/norm([nox noy]);

for n=1:N_GAUSS_POINT
    % Coordinate of the Gauss point on the reference elements
    u = GAUSS_POINT(n);
    % Shape functions for the solution
    L = [1-u 1+u]/2;
    % Shape functions for the gradient of the solution
    Lu = [-1 1]/2;
    % Jacobian to the physical element
    Js = sqrt((Lu*x_point)^2+(Lu*y_point)^2);
    % Interpolate the velocity
    Vx = L*[Vx1;Vx2];
    Vy = L*[Vy1;Vy2];
    % Pressure normal gradient
    dPdn = -1i*omega*RHO*(Vx*no(1)+Vy*no(2));
    % Add contribution to the element rhs
    Fe = Fe + GAUSS_WEIGHT(n)*( L.'*dPdn )*Js;
end

