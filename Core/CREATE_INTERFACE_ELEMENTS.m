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


function CREATE_INTERFACE_ELEMENTS(d, g, ng)
% Create a new group of elements corresponding to interfaces between
% existing elements
% d: domain number of the elements to process
% g: a vector of group numbers that should be processed
% ng: the number of the new group of interface elements

global ELEMENT N_NODE_ELEMENT N_ELEMENT
global ELEMENT_ID ELEMENT_DOMAIN

ne = GET_ELEMENT_LIST(d, 0, g);

% Table of nodes per element
e_ = ELEMENT(1:3, ne );
% Table of element edges
segment = [e_(1:2,:) e_(2:3,:) [e_(3,:);e_(1,:)]];
% Table mapping edge into original elements
seg_elem = [(1:size(e_,2)) (1:size(e_,2)) (1:size(e_,2))];
% Construct the matrix of connectivity between nodes
temp = sparse(segment(1,:),segment(2,:),seg_elem);
% Remove edges which are not double sided (remove boundary elements)
temp = temp .* ((temp~=0).*(temp.'~=0));
% Extract the connectivity between nodes
[ii,jj,ss] = find(temp);
s_ = [ii' ; jj' ; ne(ss)];
[ii,jj,ss] = find(temp.');
s_ = [s_ ; ne(ss)];
% Remove the duplicated edges
s_ = s_(:, s_(3,:)>s_(4,:) );
% Update the finite element model
ns_ = size(s_,2);
% Update the database
ELEMENT = [[ELEMENT;zeros(size(s_,1)-size(ELEMENT,1),N_ELEMENT)] s_];
N_NODE_ELEMENT = [N_NODE_ELEMENT s_(1,:)*0+2];
N_ELEMENT = N_ELEMENT + ns_;
ELEMENT_ID = [ELEMENT_ID s_(1,:)*0];
ELEMENT_DOMAIN = [ELEMENT_DOMAIN [d;0;ng]*ones(1,ns_)];
