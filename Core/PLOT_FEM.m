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


function h = PLOT_FEM(node, element, A)
% Plot a solution on a 2D finite element mesh
%
% node: a 2xN matrix where the rows corresponds to the coordinates and each
% column corresponds to a node.
%
% element: a kxM matrix where k is the number of nodes per element and M is
% the number of elements
%
% A: a vector of values. If A is of size N then it contains a value for
% each node and an interpolation of this values is used over each element.
% If A is of size M then it contains a value for each element and a
% constant value is plotted in each element.
%
% If no argument A is provided then the mesh is plotted.
%
% If requested, this function returns a handle for the patch created in the
% figure.

X = node(1,:);
Y = node(2,:);

if nargin == 2
    hp = plot(X(element), Y(element), 'b-');
elseif length(A) == size(node,2)
    hp = patch(X(element), Y(element), A(element), 'EdgeColor', 'none');
elseif length(A) == size(element,2)
    hp = patch(X(element), Y(element), A(:).', 'EdgeColor', 'none');
else
    error('Invalid input.');
end

if nargout == 1
    h = hp;
end
