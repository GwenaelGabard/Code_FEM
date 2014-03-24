function out = domain(d,x)

global ELEMENT_DOMAIN NODE ELEMENT

% domain() returns the number of domains
if nargin == 0
    out = 2;
    return
end

% domain(i) returns a description of the ith domain
if nargin == 1
    % Domain 1
    out = {{'Library\Helmholtz_1D' 'METHOD_FEM' 'L2' 'element'}};
    return
end

% domain(i,x) returns of vector of properties in the ith domain for the
% node with coordinates x
RHO0 = 1;  % Fluid mean density
C0 = 1; % Mean sound speed
if nargin == 2
    out = zeros(5,size(x,2));
    out(1,:) = RHO0;  % rho0
    out(2,:) = C0;  % c0
    out(3,:) = 0;  % U0x
    out(4,:) = 0;  % U0y
    out(5,:) = 0;  % Source term
    return
end
