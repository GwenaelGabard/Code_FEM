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


% Standard FEM assembling function
% list: a vector containing the list of the numbers of the elements that
% are part of the current domain.

fprintf('METHOD_FEM: Domain %i (%i elements)\n',ELEMENT_DOMAIN(1,list(1)),length(list));
% Calculate the maximum number of non-zero entries
nnzK = sum(N_DOF_ELEMENT(list).^2);
% Prepare temporary storage for the non-zero entries
ii = ones(nnzK,1);
jj = ii;
ss = zeros(nnzK,1);
aa = 0;
% Loop over the elements
for j=1:length(list)
    % Call the element function
    [Ke,Fe,Re,Ve] = feval(char(ELEMENT_NAME(ELEMENT_LIST(list(j)))),list(j));
    % The list of DOFs for the current element
    dof_list = DOF_ELEMENT(1:N_DOF_ELEMENT(list(j)),list(j));
    % Store the entries for the matrix
    [iie,jje,sse] = find(Ke);
    if ~isempty(iie)
        aa = aa(end) + (1:length(iie));
        ii(aa) = dof_list(iie);
        jj(aa) = dof_list(jje);
        ss(aa) = sse;
    end
    % Store the contributions to the RHS
    F(dof_list) = F(dof_list) + Fe;
    % Store the linear constraints
    Le = find(max(abs(Re),[],2));
    % Clean up any previous linear constraints
    %   R(dof_list(Le),:) = 0;
    %   R(dof_list(Le),dof_list) = Re(Le,:);
    %   V(dof_list(Le)) = Ve(Le);
    % Add up the linear constraints
    R(dof_list(Le),dof_list) = R(dof_list(Le),dof_list) + Re(Le,:);
    V(dof_list(Le)) = V(dof_list(Le)) + Ve(Le);
end
% Perform the assembling
K = K + sparse(ii,jj,ss,N_DOF,N_DOF);

% Tidy up the workspace
clear ii jj ss iie jje sse Le Fe Ke Re Ve aa
