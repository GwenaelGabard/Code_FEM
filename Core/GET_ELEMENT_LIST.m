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


function a = GET_ELEMENT_LIST(d,t,g)
% Return the list of the elements that are in domain d, and of type t
% (0: domain, 1: boundary), and in group g.
% If no group g is provided then all the elements in domain d and of type t
% are returned.

global ELEMENT_DOMAIN

if nargin==3
    a = find((ELEMENT_DOMAIN(1,:)==d)&(ELEMENT_DOMAIN(2,:)==t)&(ELEMENT_DOMAIN(3,:)==g(1)));
    if length(g)>1
        for n=2:length(g)
            a = [ a find((ELEMENT_DOMAIN(1,:)==d)&(ELEMENT_DOMAIN(2,:)==t)&(ELEMENT_DOMAIN(3,:)==g(n))) ]; %#ok<AGROW>
        end
    end
elseif nargin==2
    a = find((ELEMENT_DOMAIN(1,:)==d)&(ELEMENT_DOMAIN(2,:)==t));
end
