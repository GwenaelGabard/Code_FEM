function [Ke,Fe,Re,Ve] = potential(k)

if nargin==0
   Ke = [3 3 1 1 1 0];
   return
end

global U
global DOF_ELEMENT BOUNDARY_DATA

Ke = zeros(6);
Fe = zeros(6, 1);
Re = eye(6);
Ve = zeros(6, 1);

Ve(1:3) = BOUNDARY_DATA(1:3,k) - U(DOF_ELEMENT(1:3,k));
Re(4:6,4:6) = 0;
