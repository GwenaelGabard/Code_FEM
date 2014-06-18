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


% Construct the data describing the physics of the problem

global ELEMENT_NAME ELEMENT_LIST
global ELEMENT_DATA NODE_DATA BOUNDARY_DATA
global INTERFACE INTERFACE_NAME
global DOMAIN_METHOD DOMAIN_PATH
fprintf('BUILD_MODEL: ');

ELEMENT_NAME = [];
ELEMENT_DATA = [];
ELEMENT_LIST = zeros(1,N_ELEMENT);
NODE_DATA = sparse(1000,N_NODE);
BOUNDARY_DATA = sparse(1000,N_ELEMENT);
DOMAIN_PATH = cell(2,N_DOMAIN);
INTERFACE = cell(N_DOMAIN,N_DOMAIN);
INTERFACE_NAME = cell(N_DOMAIN,N_DOMAIN);
DOMAIN_METHOD = cell(N_DOMAIN,N_DOMAIN);
max_elem = 0;
max_node = 0;
% Loop over all the domains
for k=1:N_DOMAIN
    % Create a list of elements part of the current domain
    list = find(ELEMENT_DOMAIN(1,:)==k);
    % Call the 'domain' function to get the info about it
    temp = domain(k);
    % Store the problem name
    PROBLEM = temp{1}{1};
    % Store the method to be used for this domain
    DOMAIN_METHOD{k,k} = temp{1}{2};
    % Store the directories containing the element functions
    ZExist=0;
    for iPath = 1:length(LIBRARY_PATH)
        LocalPath = [LIBRARY_PATH{iPath} filesep PROBLEM filesep temp{1}{3} filesep];
        DOMAIN_PATH{iPath,k} = LocalPath;
        if (exist(LocalPath,'dir'))
           ZExist=1;
        end
    end
    if (ZExist == 0)
       error('Elementary folder routine could not be accessed!') 
    end
    % Loop over the elements of the domain
    for j=1:length(list)
        % If it is a domain element
        if ELEMENT_DOMAIN(2,list(j))==0
            if size(temp{1})==4
                % Store the name of the element function
                [number,ELEMENT_NAME] = ADD_TO_LIST(ELEMENT_NAME,temp{1}{4});
                ELEMENT_LIST(list(j)) = number;
            else
                % Store the name of the element function
                [number,ELEMENT_NAME] = ADD_TO_LIST(ELEMENT_NAME,temp{1}{3+ELEMENT_DOMAIN(3,list(j))});
                ELEMENT_LIST(list(j)) = number;
            end
        end
        % If it is a boundary element
        if ELEMENT_DOMAIN(2,list(j))==1
            % Call the 'boundary' function to get the description of the
            % boundary
            name = boundary(k,ELEMENT_DOMAIN(3,list(j)));
            [number,ELEMENT_NAME] = ADD_TO_LIST(ELEMENT_NAME,name);
            ELEMENT_LIST(list(j)) = number;
            name = [];
            for n=1:N_NODE_ELEMENT(list(j))
                % Store the name of the element function
                name = [name;boundary(k,ELEMENT_DOMAIN(3,list(j)),NODE(:,ELEMENT(n,list(j)))).']; %#ok<AGROW>
            end
            % Store the data associated with this boundary condition
            BOUNDARY_DATA(1:length(name),list(j)) = name; %#ok<SPRIX>
            max_elem = max(max_elem,length(name));
        end
    end
    % If there is coupling between domains
    for j=2:size(temp,1)
        DOMAIN_METHOD{k,temp{j}{1}} = temp{j}{2};
        INTERFACE_NAME{k,temp{j}{1}} = temp{j}{3};
        INTERFACE{k,temp{j}{1}} = temp{j}{4};
    end
    % Call the 'domain' function to get the data associated with the nodes
    list = find(NODE_DOMAIN==k);
    temp = domain(k,NODE(:,list));
    NODE_DATA(1:size(temp,1),list) = temp; %#ok<SPRIX>
    max_node = max(max_node,size(temp,1));
end

NODE_DATA = full(NODE_DATA(1:max_node,:));
BOUNDARY_DATA = full(BOUNDARY_DATA(1:max_elem,:));

fprintf('%i domains, %i nodes, %i elements\n',N_DOMAIN,N_NODE,N_ELEMENT);

% Remove temporary variables
clear k j n temp list number name ii jj ss
clear max_node max_elem
clear public private