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


function out = domain(d,x)

global ELEMENT_DOMAIN NODE ELEMENT

% domain() returns the number of domains
if nargin == 0
    out = 2;
    return
end

% domain(i) returns a description of the ith domain
if nargin == 1
    % Domain 1
    out = {{'Library/Helmholtz_2D' 'METHOD_FEM' 'T3' 'element'}};
    return
end

% domain(i,x) returns of vector of properties in the ith domain for the
% node with coordinates x
RHO0 = 1;  % Fluid mean density
C0 = 1; % Mean sound speed
if nargin == 2
    out = zeros(5,size(x,2));
    out(1,:) = RHO0;  % rho0
    out(2,:) = C0;  % c0
    out(3,:) = 0;  % U0x
    out(4,:) = 0;  % U0y
    out(5,:) = 0;  % Source term
    return
end
