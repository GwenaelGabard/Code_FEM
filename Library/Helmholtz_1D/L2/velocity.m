function [Ke,Fe,Re,Ve] = velocity(k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Convected Helmholtz Equation
% Notation: exp(-i \omega t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Point element
% 1 node
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
    Ke = [1 1 1 0];
    return
end

global omega
global ELEMENT NODE
global  NODE_DATA

RHO = NODE_DATA(1,ELEMENT(1,k));
C0 = NODE_DATA(2,ELEMENT(1,k));

Ke = 0;
Fe = 0;
Re = 0;
Ve = 0;

% RHS, piston
Fe = +1i*(omega/C0);
