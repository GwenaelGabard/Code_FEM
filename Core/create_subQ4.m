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


function [u,v,c] = create_subQ4(N)
% Create a set of nodes and small triangular elements covering the
% reference quadrilateral element (with 0<u<1 and 0<v<1).
% This can be used to construct a refined interpolation of the solution on
% an element.
% 
% Example:
% [u,v,c] = create_subQ4(N);
% plot(u,v,'ok',u(c).',v(c).','r-');
% 
% N: number of elements created along the side of the reference element.
% u, v: vectors of coordinates of the nodes.
% c: matrix containing the 3 nodes for each triangle, with each row
% corresponding to a triangle.


% Create the nodes
u = linspace(0,1,N+1);
v = u;
[u,v] = meshgrid(u,v);
u = u.';
v = v.';
u = u(:);
v = v(:);

% Create the elements
c = zeros(N^2,3);
k = 1;
for nv=1:N
    for nu=1:N
        c(k,:) = nu + (nv-1)*(N+1) + [0 1 N+1];
        k = k+1;
        c(k,:) = nu + (nv-1)*(N+1) + [1 N+2 N+1];
        k = k+1;
    end
end
