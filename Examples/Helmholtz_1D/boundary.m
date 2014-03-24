function [out] = boundary(d,bs,x)

% boundary(d, bs) returns the type of boundary condition for the boundary
% bs of domain d
if nargin == 2
    list = {'velocity','impedance'};
    out = list{bs};
    return
end

global omega

% boundary(d, bs, x) returns the properties associated to the point x on
% the boundary bs of domain d
if nargin == 3
    switch bs
        case 1
            out = 0;
        case 2
            out = 1i*omega;
    end
end
