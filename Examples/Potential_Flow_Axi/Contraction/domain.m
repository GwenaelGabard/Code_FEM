function [out] = domain(~,x)
% domain() returns the number of domains
% domain(i) returns a description of the ith domain
% domain(i,x) returns of vector of properties in the ith domain for the
% node with coordinates x

if nargin == 0
    out = 1;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 1
    out = {{'Library/Potential_Flow_Axi' 'METHOD_FEM' 'T6' 'element'}};
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 2
    out = zeros(1,size(x,2));
    return
end
