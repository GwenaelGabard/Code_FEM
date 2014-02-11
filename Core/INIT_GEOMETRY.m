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


% Declare and initialise the variables needed to describe the geometry of
% the problem

clear N_DIM N_NODE NODE
clear N_ELEMENT N_NODE_ELEMENT ELEMENT
clear N_DOMAIN ELEMENT_DOMAIN NODE_DOMAIN

global N_DIM N_NODE NODE %#ok<NUSED>
global N_ELEMENT N_NODE_ELEMENT ELEMENT
global N_DOMAIN ELEMENT_DOMAIN NODE_DOMAIN
global ELEMENT_ID NODE_ID

N_NODE = 0;
N_ELEMENT = 0;
N_DOMAIN = 0;
NODE = [];
N_NODE_ELEMENT = [];
ELEMENT = [];
ELEMENT_DOMAIN = [];
NODE_DOMAIN = [];
ELEMENT_ID = [];
NODE_ID = [];
