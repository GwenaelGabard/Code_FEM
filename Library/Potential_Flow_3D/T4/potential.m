function [Ke,Fe,Re,Ve] = potential(k)

if nargin==0
   Ke = [3 3 1 1 1 0];
   return
end

global U
global DOF_ELEMENT BOUNDARY_DATA

Ke = zeros(4);
Fe = zeros(4,1);
Re = eye(4);
Ve = zeros(4,1);

Ve(1:3) = BOUNDARY_DATA(1:3,k) - U(DOF_ELEMENT(1:3,k));
Re(4,4) = 0;