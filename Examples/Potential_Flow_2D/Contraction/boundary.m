function [out] = boundary(~,bs,~)
% boundary(d, bs) returns the type of boundary condition for the boundary
% bs of domain d
% boundary(d, bs, x) returns the properties associated to the point x on
% the boundary bs of domain d

% Boundaries 1, 2, 3: lower wall (straight, curved, straight)
% Boundary 4: duct termination on the right
% Boundaries 5, 6, 7: upper wall (straight, curved, straight)
% Boundary 8: duct termination on the left

if nargin == 2
    list = {'velocity','velocity','velocity','velocity','velocity','velocity','velocity','potential'};
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
            out = 0;
        case 3
            out = 0;
        case 4
            out = Flow.outlet_velocity;
        case 5
            out = 0;
        case 6
            out = 0;
        case 7
            out = 0;
        case 8
            out = 0;
    end
end
