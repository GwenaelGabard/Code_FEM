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


function [Ke,Fe,Re,Ve] = dirichlet(k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Convected Helmholtz Equation
% Notation: exp(-i \omega t)
% Boundary condition: dirichlet
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
% 1: solution at node 1
% 2: solution at node 2

if nargin==0
    % Return the properties of the element
    % Ke should contain:
    % the priority of the element
    % the number of nodes
    % the number of DOFs for node 1
    % the number of DOFs for node 2
    % ...
    % the number of internal DOFs of the element
    Ke = [4 2 1 1 0];
    return
end

global BOUNDARY_DATA

Ke = zeros(2,2);
Fe = zeros(2,1);
Re = eye(2,2);
Ve = BOUNDARY_DATA(1:2,k);
