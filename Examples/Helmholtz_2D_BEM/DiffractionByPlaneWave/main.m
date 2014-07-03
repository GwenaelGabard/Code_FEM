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

% To find the solver
addpath('../../../../Code_FEM/Core');
MANAGE_PATH

global omega typegeo type Order facteurdensite
global alpha beta
global LENGTH_ELEMENT NGPWL
global N_GAUSS_POINT GAUSS_POINT GAUSS_WEIGHT SHAPEFUNC DSHAPEFUNC

% Diffraction of a circle by an incident Plane Wave

%%%%%%%%%% Plane wave Excitation  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Type of geometry
typegeo = 'circle';
Perimeter = 2*pi;
% Plane wave amplitude
Ainc = 1;
% Frequency
ka = 6;
omega = ka*340;
% Dof density
ndof = 30;
% Interpolation order
Order = 2;
% Number of BEM elements
Nelement = round(ndof*Perimeter*(ka/(2*pi))/Order);
% Number of Gauss Points per wavelength
NGPWL = 10;
% Density factor (For the auto-influence computation)
facteurdensite = 1;
% Burton and Miller
if strcmp(type,'open')
    alpha = 1; beta = 0;
else
    alpha= 1i/ka; beta = 1;
end
% Mesh creation
INIT_GEOMETRY
[N_DIM,N_NODE,NODE,N_ELEMENT,N_NODE_ELEMENT,ELEMENT,N_DOMAIN...
    ,ELEMENT_DOMAIN,NODE_DOMAIN,ELEMENT_ID,NODE_ID] =  feval('geobem',Nelement,Order);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the element length
LENGTH_ELEMENT = ElementLength(Order);
ndofpwl = (N_NODE)/((sum(LENGTH_ELEMENT))/(2*pi/ka));
% Precompute the shape functions at the Gauss points
N_GAUSS_POINT = max(max(ceil(LENGTH_ELEMENT*ka/(2*pi)*NGPWL),4));
[GAUSS_POINT,GAUSS_WEIGHT] = gauleg(-1,1,N_GAUSS_POINT);
SHAPEFUNC = zeros(Order+1,N_GAUSS_POINT);
DSHAPEFUNC = zeros(Order+1,N_GAUSS_POINT);
for n=1:N_GAUSS_POINT
    [Nu,dNu]=shape_function(GAUSS_POINT(n),Order);
    SHAPEFUNC(:,n)=Nu;
    DSHAPEFUNC(:,n)=dNu;
end
fprintf('\n Numerical data: \n ka = %i, %1.2f dof per wavelength, Order of interp. = %i \n \n',ka,ndofpwl,Order);

BUILD_MODEL
BUILD_DOF
BUILD_SYSTEM
SOLVE_SYSTEM

[Uanal] = CirclePlanewave_AnalyticalSolution(NODE,1);
ErL2 = L2Error(U,Uanal);
fprintf('L2 error: %g pc \n', ErL2);

% plot the result
PlotResults