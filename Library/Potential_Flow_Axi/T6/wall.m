function [Ke,Fe,Re,Ve] = wall(~)

if nargin==0
    Ke = [4 3 1 1 1 0];
    return
end

Ke = zeros(6);
Fe = zeros(6, 1);
Re = zeros(6);
Ve = zeros(6, 1);
