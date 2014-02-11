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


function out=POST_mesh(k)
% k is the number of the domain to plot

global ELEMENT NODE ELEMENT_DOMAIN

% Return the name of the button shown in the GUI
if nargin==0
   out = 'Mesh';
   return
end

% List of the elements in the domain
elem_list = find((ELEMENT_DOMAIN(1,:)==k)&(ELEMENT_DOMAIN(2,:)==0));
% List of the border elements in the domain
border_list = find((ELEMENT_DOMAIN(1,:)==k)&(ELEMENT_DOMAIN(2,:)==1));

% Nodes coordinates
X = NODE(1,:);
Y = NODE(2,:);
X = full(X(ELEMENT(1:3,elem_list)));
Y = full(Y(ELEMENT(1:3,elem_list)));
% Number of boundaries
n_border = max(ELEMENT_DOMAIN(3,border_list));

figure
flagmap = [0 0 1;0 1 0;1 0 0;1 1 0;0 1 1;1 0 1;0 0 0;1 1 1];
hold on
patch(X,Y,ones(size(elem_list))*'b','EdgeColor',flagmap(1,:),'FaceColor','none');
for j=1:n_border
   X = NODE(1,:);
   Y = NODE(2,:);
   border_list = find((ELEMENT_DOMAIN(1,:)==k)&(ELEMENT_DOMAIN(2,:)==1)&(ELEMENT_DOMAIN(3,:)==j));
   X = full(X(ELEMENT(1:2,border_list)));
   Y = full(Y(ELEMENT(1:2,border_list)));
   patch(X,Y,'w','EdgeColor',flagmap(j+1,:));
end
hold off
axis equal;axis tight;
