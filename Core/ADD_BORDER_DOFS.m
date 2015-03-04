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


function ADD_BORDER_DOFS(nd, border_list, group_list)
% Modify the degrees of freedom of a boundary element to include the
% degrees of freedom of the adjacent domain element.
% nd: domain number
% border_list: list of border numbers that should be modified
% group_list: list of group numbers for the domain elements to consider

global ELEMENT N_NODE_ELEMENT
global N_DOF_ELEMENT DOF_ELEMENT

DOF_ELEMENT = sparse(DOF_ELEMENT);

elem_list = GET_ELEMENT_LIST(nd,1,border_list);

if nargin==3
    dom_list = GET_ELEMENT_LIST(nd,0,group_list);
else
    dom_list = GET_ELEMENT_LIST(nd,0);
end

for n=1:length(elem_list)
    table = zeros(1,length(dom_list));
    for m=1:N_NODE_ELEMENT(elem_list(n))
        table = table + sum(ELEMENT(:,dom_list)==ELEMENT(m,elem_list(n)),1);
    end
    % Adjacent element number
    num = find(table==N_NODE_ELEMENT(elem_list(n)));
    % Problem detection
    if length(num)~=1
        error('A boundary element does not have a neighbour of has multiple ones !');
    end
    % Loop over the dofs of the adjacent element
    for m=1:N_DOF_ELEMENT(dom_list(num))
        % If the dof is already in the list, do not add it
        if ~any(DOF_ELEMENT(1:N_DOF_ELEMENT(elem_list(n)),elem_list(n))==DOF_ELEMENT(m,dom_list(num)))
            N_DOF_ELEMENT(elem_list(n)) = N_DOF_ELEMENT(elem_list(n)) + 1;
            DOF_ELEMENT(N_DOF_ELEMENT(elem_list(n)),elem_list(n)) = DOF_ELEMENT(m,dom_list(num));
        end
    end
    % Loop over the nodes of the neighbouring element
    k = N_NODE_ELEMENT(elem_list(n));
    for m=1:N_NODE_ELEMENT(dom_list(num))
        % If the node is already in the list, do not add it
        if ~any(ELEMENT(1:k,elem_list(n))==ELEMENT(m,dom_list(num)))
            k = k+1;
            ELEMENT(k,elem_list(n)) = ELEMENT(m,dom_list(num));
        end
        N_NODE_ELEMENT(elem_list(n)) = k;
        ELEMENT(k+1,elem_list(n)) = dom_list(num);
    end
end

% Back to standard coding
DOF_ELEMENT = full(DOF_ELEMENT);

% Commentaires
% This function was tested for the T3, L2, L3 and T4 elements
