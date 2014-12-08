function [node,edge,element] = convert_T3_T6(node,edge,element,geom)

side = edge;
elem = element;

segment = [elem(1:2,:) elem(2:3,:) [elem(3,:);elem(1,:)]];
seg_elem = [(1:size(elem,2)) (1:size(elem,2)) (1:size(elem,2))];
temp = sparse(segment(1,:),segment(2,:),seg_elem);
[ii,jj,ss] = find(temp .* ((temp~=0).*(temp.'~=0)));
edge = [ii' ; jj' ; ss'];
[~,~,ss] = find(temp.' .* ((temp~=0).*(temp.'~=0)));
edge = [edge ; ss'];
edge = edge(:,edge(3,:)>edge(4,:));
[ii,jj,ss] = find(temp .* not((temp~=0).*(temp.'~=0)));
edge = [edge [ii' ; jj' ; ss' ; 0*ss' ]];

X = node(1,:);
Y = node(2,:);

new_node = [ X(edge(1,:))+X(edge(2,:)) ; Y(edge(1,:))+Y(edge(2,:)) ]/2;
n_new_node = size(new_node,2);
n_node = size(node,2);
node = [ node new_node ];

% Add the new nodes to the definitions of the elements
element = zeros(6,size(elem,2));
element([1 3 5],:) = elem(1:3,:);
for n=1:n_new_node
    n1 = edge(1,n);
    n2 = edge(2,n);

    e = edge(3,n);
    if ((elem(1,e)==n1)&&(elem(2,e)==n2))||((elem(1,e)==n2)&&(elem(2,e)==n1))
        element(2,e) = n + n_node;
    end
    if ((elem(2,e)==n1)&&(elem(3,e)==n2))||((elem(2,e)==n2)&&(elem(3,e)==n1))
        element(4,e) = n + n_node;
    end
    if ((elem(3,e)==n1)&&(elem(1,e)==n2))||((elem(3,e)==n2)&&(elem(1,e)==n1))
        element(6,e) = n + n_node;
    end

    e = edge(4,n);
    if e~=0
        if ((elem(1,e)==n1)&&(elem(2,e)==n2))||((elem(1,e)==n2)&&(elem(2,e)==n1))
            element(2,e) = n + n_node;
        end
        if ((elem(2,e)==n1)&&(elem(3,e)==n2))||((elem(2,e)==n2)&&(elem(3,e)==n1))
            element(4,e) = n + n_node;
        end
        if ((elem(3,e)==n1)&&(elem(1,e)==n2))||((elem(3,e)==n2)&&(elem(1,e)==n1))
            element(6,e) = n + n_node;
        end
    end
end
element = [ element ; elem(end,:) ];

ns = find(edge(4,:)==0);
edge = edge(1:4,ns);
side = [ side ; zeros(2,size(side,2)) ];
side(end-2,:) = side(5,:);
side(8:9,:) = side(3:4,:);
for n=1:length(ns)
    n1 = edge(1,n);
    n2 = edge(2,n);

    for m=1:size(side,2)
        if ((side(1,m)==n1)&&(side(2,m)==n2))||((side(1,m)==n2)&&(side(2,m)==n1))
            side(4,m) = ns(n) + n_node;
            side(5,m) = edge(3,n);
        end
    end
end
edge = side([1 4 2 3 5 7 8 9],:);

% If provided use the geometry to move the nodes on the boundary exactly on
% the surface.
if nargin==4
    for n=1:size(edge,2)
        nn = edge(2,n);
        [x,y] = feval(geom, edge(6,n), mean(side(8:9,n)) );
        node(1,nn) = x;
        node(2,nn) = y;
    end
end

% Node renumbering to reduce the bandwidth
[p,inv_p] = reorder_node(element(1:6,:),6);
node = node(:,p);
element(1:6,:) = inv_p(element(1:6,:));
temp = edge(1:3,:);
edge(1:3,:) = inv_p(temp);
