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

function [Ke,Fe,Re,Ve] = empty(k)


if nargin==0
    % Return the properties of the element
    % Ke should contain:
    % the priority of the element
    % the number of nodes
    % the number of DOFs for node 1
    % the number of DOFs for node 2
    % ...
    % the number of internal DOFs of the element
   Ke = [1 2 1 1 0];
   return
end

Ke = zeros(2,2);
Fe = zeros(2,1);
Re = zeros(2,2);
Ve = zeros(2,1);
