function [out] = boundary(~,bs,~)
% boundary(d, bs) returns the type of boundary condition for the boundary
% bs of domain d
% boundary(d, bs, x) returns the properties associated to the point x on
% the boundary bs of domain d

% Boundaries 1: cylindrical axis
% Boundary 2: duct termination on the right
% Boundaries 3, 4, 5: upper wall (straight, curved, straight)
% Boundary 6: duct termination on the left

if nargin == 2
    list = {'velocity','velocity','velocity','velocity','velocity','potential'};
    out = list{bs};
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Flow
if nargin == 3
    switch bs
        case 1
            out = 0;
        case 2
            out = Flow.outlet_velocity;
        case 3
            out = 0;
        case 4
            out = 0;
        case 5
            out = 0;
        case 6
            out = 0;
    end
end
