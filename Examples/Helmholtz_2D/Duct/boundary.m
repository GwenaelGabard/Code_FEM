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


function [out] = boundary(d,bs,x)

% boundary(d, bs) returns the type of boundary condition for the boundary
% bs of domain d
if nargin == 2
    list = {'wall', 'impedance', 'wall', 'wall'};
    out = list{bs};
    return
end

% boundary(d, bs, x) returns the properties associated to the point x on
% the boundary bs of domain d
if nargin == 3
    switch bs
        case 1
            out = 0;
        case 2
            out = 1;
        case 3
            out = 0;
        case 4
            out = 1;
    end
end
