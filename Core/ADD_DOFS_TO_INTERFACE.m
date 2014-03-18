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


function ADD_DOFS_TO_INTERFACE(d, t, g)
% Modify interface elements so that they include the degrees of freedom of
% the two elements on either side.
% d: domain number
% t: type of elements
% g: list of group numbers of the interface elements to modify

global ELEMENT N_NODE_ELEMENT
global DOF_ELEMENT N_DOF_ELEMENT

list = GET_ELEMENT_LIST(d,t,g);

for e=list
    % Number of the elements on either side
    e1 = ELEMENT(N_NODE_ELEMENT(e)+1,e);
    e2 = ELEMENT(N_NODE_ELEMENT(e)+2,e);
    % Number of degrees of freedom on either side
    n1 = N_DOF_ELEMENT(e1);
    n2 = N_DOF_ELEMENT(e2);
    % Update the interface element
    N_DOF_ELEMENT(e) = n1 + n2;
    DOF_ELEMENT(1:N_DOF_ELEMENT(e),e) = [DOF_ELEMENT(1:n1,e1);DOF_ELEMENT(1:n2,e2)];
end
