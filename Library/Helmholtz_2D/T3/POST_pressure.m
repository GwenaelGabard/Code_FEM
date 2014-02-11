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


function out = POST_pressure(k)
% k is the number of the domain to plot

global U DOF_ELEMENT NODE ELEMENT_DOMAIN ELEMENT

% Return the name of the button shown in the GUI
if nargin==0
   out = 'Pressure';
   return
end

% List of the element in the domain
elem_list = find((ELEMENT_DOMAIN(1,:)==k)&(ELEMENT_DOMAIN(2,:)==0));

% Extract the field to plot
DATA = full(U(DOF_ELEMENT(1:3,elem_list)));
DATA = real(DATA);

% Coordinates of the nodes
X = NODE(1,:);
Y = NODE(2,:);
X = full(X(ELEMENT(1:3,elem_list)));
Y = full(Y(ELEMENT(1:3,elem_list)));

% Draw the elements
opengl('NeverSelect');
figure
patch(X,Y,DATA,'Edgecolor','none','Erasemode','normal');
MIN = min(min(DATA));
MAX = max(max(DATA));
caxis([MIN MAX]);
colorbar;
axis equal; axis tight; box on;
xlabel('X');
ylabel('Y');
title('Pressure');
