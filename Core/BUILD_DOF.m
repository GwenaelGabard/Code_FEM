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


% Construct lists of degrees of freedom
global N_DOF_ELEMENT N_IDOF_ELEMENT N_DOF
global DOF_NODE N_DOF_NODE DOF_ELEMENT
fprintf('BUILD_DOF: ');

% N_IDOF_ELEMENT(i) : number of internal degrees of freedom for the ith
% element
N_IDOF_ELEMENT = zeros(1,N_ELEMENT);
% N_DOF_ELEMENT(i) : number of degrees of freedom for the ith element
N_DOF_ELEMENT = zeros(1,N_ELEMENT);
% The level of priority of the elements
priority = zeros(1,N_ELEMENT);
% Go through each domain at a time
for k=1:N_DOMAIN
    % Create a list of elements for this domain
    list = find(ELEMENT_DOMAIN(1,:)==k);
    % Get access to the element functions for this domain
    if (exist(DOMAIN_PATH{1,k},'dir'))
       addpath(DOMAIN_PATH{1,k});
    end
    if (exist(DOMAIN_PATH{2,k},'dir'))
       addpath(DOMAIN_PATH{2,k});
    end
    for j=1:length(list)
        % Call the element function
        temp = feval(char(ELEMENT_NAME(ELEMENT_LIST(list(j)))));
        % Check that the number of nodes given by the element function is
        % consistent with what is given in the geometry
        if N_NODE_ELEMENT(list(j))~=temp(2)
            error('Incompatibility between the element geometrical description and the element assembly function.');
        end
        % Number of internal degrees of freedom
        N_IDOF_ELEMENT(list(j)) = temp(end);
        % Number of degrees of freedom 
        N_DOF_ELEMENT(list(j)) = sum(temp(3:end));
        % The level of priority of this element compared to the others
        priority(list(j)) = temp(1);
    end
    if (exist(DOMAIN_PATH{1,k},'dir'))
       rmpath(DOMAIN_PATH{1,k});
    end
    if (exist(DOMAIN_PATH{2,k},'dir'))
       rmpath(DOMAIN_PATH{2,k});
    end
end

% Associate degrees of freedom to nodes and elements
% N_DOF_NODE(i) is the number of degrees of freedom for the ith node
N_DOF_NODE = zeros(1,N_NODE);
% N_DOF is the total number of degrees of freedom
N_DOF = 0;
% DOF_NODE(i,j) is the number of the ith degree of freedom of the jth node
DOF_NODE = zeros(max(N_DOF_NODE),N_NODE);
% DOF_ELEMENT(i,j) is the number of the ith degree of freedom of the jth
% element
DOF_ELEMENT = zeros(max(N_DOF_ELEMENT),N_ELEMENT);
% Work each domain at a time
for k=1:N_DOMAIN
    % Create a list of elements for this domain
    list = find(ELEMENT_DOMAIN(1,:)==k);
    % Create a list of nodes for this domain
    node_list = find(NODE_DOMAIN==k);
    % Get access to the element functions for this domain
    if (exist(DOMAIN_PATH{1,k},'dir'))
       addpath(DOMAIN_PATH{1,k});
    end
    if (exist(DOMAIN_PATH{2,k},'dir'))
       addpath(DOMAIN_PATH{2,k});
    end
    % Loop over all the elements of this domain
    for j=1:length(list)
        % Call element function to know the number of degrees of freedom
        temp = feval(char(ELEMENT_NAME(ELEMENT_LIST(list(j)))));
        % Store the number of degrees of freedom for the nodes
        N_DOF_NODE(ELEMENT(1:N_NODE_ELEMENT(list(j)),list(j))) = temp(3:end-1);
    end
    % Loop over the nodes of the domain
    for j=1:length(node_list)
        % Number the degrees of freedom and store these numbers
        DOF_NODE(1:N_DOF_NODE(node_list(j)),node_list(j)) = ( (N_DOF+1):(N_DOF+N_DOF_NODE(node_list(j))) )';
        N_DOF = N_DOF+N_DOF_NODE(node_list(j));
    end
    % Loop over all the elements of this domain
    for j=1:length(list)
        % Associate the degrees of freedom to the elements
        temp = DOF_NODE(:,ELEMENT(1:N_NODE_ELEMENT(list(j)),list(j)));
        temp = temp(:);
        DOF_ELEMENT(1:N_DOF_ELEMENT(list(j))-N_IDOF_ELEMENT(list(j)),list(j)) = temp(find(temp));
        DOF_ELEMENT(N_DOF_ELEMENT(list(j))-N_IDOF_ELEMENT(list(j))+1:N_DOF_ELEMENT(list(j)),list(j)) = ( (N_DOF+1):(N_DOF+N_IDOF_ELEMENT(list(j))) )';
        N_DOF = N_DOF + N_IDOF_ELEMENT(list(j));
    end
    if (exist(DOMAIN_PATH{1,k},'dir'))
       rmpath(DOMAIN_PATH{1,k});
    end
    if (exist(DOMAIN_PATH{2,k},'dir'))
       rmpath(DOMAIN_PATH{2,k});
    end
end
fprintf('%i DOFs\n',N_DOF);

% Remove temporary variables
clear order temp priority k j list node_list
