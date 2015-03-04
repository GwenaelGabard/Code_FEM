function [Ke,Fe,Re,Ve] = wall(~)

if nargin==0
    Ke = [4 3 1 1 1 0];
    return
end

Ke = zeros(4);
Fe = zeros(4, 1);
Re = zeros(4);
Ve = zeros(4, 1);
