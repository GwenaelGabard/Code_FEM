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


% This script drives the assembly of the global system

% Initialise an empty sparse matrix for the global matrix
K = sparse(N_DOF,N_DOF);
% Initialise an empty sparse vector for the right-hand side
F = sparse(N_DOF,1);
% Initialise an empty sparse matrix for the linear constraints
R = sparse(N_DOF,N_DOF);
% Initialise an empty sparse vector for the linear constraints
V = sparse(N_DOF,1);
% Loop over the domains
for k=1:N_DOMAIN
    % Create a list of the elements in the current domain
    list = find(ELEMENT_DOMAIN(1,:)==k); %#ok<NASGU>
    % Add the path to the element functions for this domain
    addpath(char(DOMAIN_PATH{k}));
    % Launch the assembly script for this domain
    eval(char(DOMAIN_METHOD{k,k}));
    % Go through possible couplings with other domains
    for j=1:N_DOMAIN
        % Check if the domain k and j are coupled
        if ~isempty(INTERFACE{k,j})
            % Create the list of interactions between elements of the two
            % domains
            list = INTERFACE{k,j}; %#ok<NASGU>
            % Launch the script to calculate the interaction
            eval(char(DOMAIN_METHOD{k,j}));
        end
    end
    % Remove the path to the element functions for this domain
    rmpath(char(DOMAIN_PATH{k}));
end

% Remove temporary variables
clear list k j
