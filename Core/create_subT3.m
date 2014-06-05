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


function [u,v,c] = create_subT3(N)
% Create a set of nodes and small triangular elements covering the
% reference triangular element (with 0<u<1 and 0<v<1).
% This can be used to construct a refined interpolation of the solution on
% an element.
% 
% Example:
% [u,v,c] = create_subT3(N);
% plot(u,v,'ok',u(c).',v(c).','r-');
% 
% N: number of elements created along the side of the reference element.
% u, v: vectors of coordinates of the nodes.
% c: matrix containing the 3 nodes for each triangle, with each row
% corresponding to a triangle.


% Create the nodes
u = zeros(1,sum(N+1:-1:1));
v = u;
k = 1;
for nv=0:N
    for nu=0:(N-nv)
        u(k) = nu/N;
        v(k) = nv/N;
        k = k+1;
    end
end


% Create the elements
if N==1
    c = [1 2 3];
    return
end
c = zeros(sum(2*(N:-1:1)-1),3);
k = 1;
for nv=1:N
    for nu=1:(N-nv+1)
        if (nv==1)
            if (nu==1)
                c(1,:) = [1 2 N+2];
                c(2,:) = [2 N+3 N+2];
                k = 3;
            else
                c(k,:) = c(k-2,:)+1;
                k = k+1;
                if nu<(N-nv+1)
                    c(k,:) = c(k-2,:)+1;
                    k = k+1;
                end
            end
        else
            if (nu==1)
                c(k,:) = c(k-(N-nv+2)*2+1,1) + N+3-nv + [ 0 1 N-nv+2];
                k = k+1;
                if nu<(N-nv+1)
                    c(k,:) = c(k-(N-nv+2)*2+1,1) + N+3-nv + [ 0 N-nv+2 N-nv+1];
                    k = k+1;
                end
            else
                c(k,:) = c(k-2,:)+1;
                k = k+1;
                if nu<(N-nv+1)
                    c(k,:) = c(k-2,:)+1;
                    k = k+1;
                end
            end
        end
    end
end
