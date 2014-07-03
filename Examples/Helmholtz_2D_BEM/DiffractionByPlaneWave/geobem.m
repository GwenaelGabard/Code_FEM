function [n_dim,n_node,node,n_element,n_node_element,element,n_domain,element_domain,node_domain,element_ID,node_ID] = geobem(n_element,Order)

global typegeo type

% Generates a BEM 2D geometry
% Order designates the interpolation order (Lagrange)

% Number of central nodes
Nc = Order-1;

% Number of nodes
n_node = n_element*(Nc+1);

% Angle parametrization
t = linspace(0,2*pi,n_node+1); t = t(1:end-1);

switch typegeo
    % CLOSED BOUNDARIES
    case 'circle'
        type = 'closed';
        X = cos(t);
        Y = sin(t);
    case 'ellipse'
        type = 'closed';
        a = 2;  b = 1;
        X = a*cos(t); Y = b*sin(t);
    case 'boomerang'
        type = 'closed';
        X = (cos(t)+10*cos(2*t)-10)/10;
        Y = sin(t);
    case 'star'
        type = 'closed';
        X = cos(t).*(1+0.5*cos(10*t));
        Y = sin(t).*(1+0.5*cos(10*t));
    case 'square'
        type = 'closed';
        if (n_element/4)~=round(n_element/4)%
            fprintf('Multiple of 4 is preffered');
            n_element = round(n_element/4)*4;
            n_node = n_element*(Nc+1);
        end
        L = 2;
        nnc = n_node/4;
        X(1:nnc) = (L/2)-L/(nnc)*(0:nnc-1); Y(1:nnc) = 1;
        X(nnc+1:2*nnc) = -1; Y(nnc+1:2*nnc) = (L/2)-L/(nnc)*(0:nnc-1);
        X(2*nnc+1:3*nnc) = -(L/2)+L/(nnc)*(0:nnc-1); Y(2*nnc+1:3*nnc) = -1;
        X(3*nnc+1:4*nnc) = 1; Y(3*nnc+1:4*nnc) = -(L/2)+L/(nnc)*(0:nnc-1);
    
    % OPEN BOUNDARIES
    case 'circlearc'
        type = 'open';
        n_node = n_node+1;
        th = linspace(pi/2,3*pi/2,n_node);
        X = cos(th);
        Y = sin(th);
        
    case 'plate'
        type = 'open';
        n_node = n_node+1;
        L = 2;
        Y(1:n_node) = (L/(n_node-1)*(0:n_node-1))-L/2; X(1:n_node) = 0;
        
    case 'spiral'
        type = 'open';
        n_node = n_node+1;
        t = linspace(0,2*pi,n_node);
        X = (t.*cos(2*t))/(2*pi);
        Y = (t.*sin(2*t))/(2*pi);
        
    otherwise
        error('Geometry not defined');
end

node(1,:) = X;
node(2,:) = Y;
n_node_element = Nc+2;

% Connectivity creation
element(1,:) = 1:Nc+1:n_node-Nc;
for n=2:n_node_element
    element(n,:) = element(n-1,:)+1;
end
if (strcmp(type,'open')==0)
    element(n_node_element(end),end) = element(1,1);
end

n_dim = 2;
n_domain = 1;
element_domain = [ones(1,n_element); zeros(1,n_element); ones(1,n_element)] ; % QUE DES elementS DE DOMAINE
node_domain = ones(1,n_node);
element_ID = 1:n_element;
node_ID = 1:n_node;
n_node_element = n_node_element+zeros(1,n_element);

