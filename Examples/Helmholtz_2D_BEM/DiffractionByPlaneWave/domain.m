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

function [out] = domain(d,x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0
    out = 1;
    return
end

if nargin == 1
    out = {{'Library/Helmholtz_2D_BEM' 'METHOD_BEM' 'L2' 'DBEM_Wall_Wall'}};
    return
end

RHO0 = 1.2;
C0 = 340;

if nargin == 2
    out = RHO0*ones(1,size(x,2));
    out = [out ; C0*ones(1,size(x,2))];
    
    return
end
