function [Ke,Fe,Re,Ve] = element(k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Convected Helmholtz Equation
% Notation: exp(-i \omega t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear element
% 2 nodes
% Data in ELEMENT_DATA(:,k):
% None
% Data in NODE_DATA(:,k):
% rho0, c0, u0, v0, s

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

global omega
global ELEMENT NODE
global  NODE_DATA

RHO = NODE_DATA(1,ELEMENT(1,k));
C0 = NODE_DATA(2,ELEMENT(1,k));

x_node = NODE(1,ELEMENT(1:2,k));

Ke = zeros(2,2);
Fe = zeros(2,1);
Re = zeros(2,2);
Ve = zeros(2,1);

N_GAUSS_POINT = 2;
GAUSS_POINT = [-1/sqrt(3) 1/sqrt(3)];
GAUSS_WEIGHT = [1 1];

for n=1:N_GAUSS_POINT
   u = GAUSS_POINT(n);
   Le = x_node(2)-x_node(1);
   N = [(1-u)/2 (1+u)/2];
   GRAD = 2/Le*[-1/2 1/2];
   Ke = Ke + GAUSS_WEIGHT(n)*( GRAD.'*GRAD - (omega/C0)^2*N.'*N )*Le/2;  
end


