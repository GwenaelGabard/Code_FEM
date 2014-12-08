function [p,inv_p] = reorder_element(element,nn)
%
% Calculate an optimized numbering of the elements that minimizes the
% bandwidth of the element graph.
% Usage:
%   [p,inv_p] = reorder_element(element,nn)
% element is a matrix where each column corresponds to an element.
% nn is the number of nodes per element.
% p is a permutation vector.
% p_inv is the inverse permutation vector.

% Compute the connectivity between elements
segment = [ ];
seg_elem = 1:size(element,2);
for n=1:nn-1
    segment = [ segment element(n:n+1,:) ]; %#ok<AGROW>
    seg_elem = [ seg_elem (1:size(element,2)) ]; %#ok<AGROW>
end
segment = [ segment [ element(nn,:) ; element(1,:) ] ];
temp = sparse(segment(1,:),segment(2,:),seg_elem);
[ii,jj,ss] = find(temp .* ((temp~=0).*(temp.'~=0)));
edge = [ii' ; jj' ; ss'];
[~,~,ss] = find(temp.' .* ((temp~=0).*(temp.'~=0)));
edge = [edge ; ss'];
edge = edge(:,edge(3,:)>edge(4,:));

% Build a spare matrix with 1 on the component a_ij when the elements i and
% j are neighbours.
K = sparse(edge(3,:),edge(4,:),ones(1,size(edge,2)),size(element,2),size(element,2));
% Use the symetric reversed Cuthill-McKee method to reduce the bandwidth
p = symrcm(K);
% Get the inverse permutation
[~,inv_p] = sort(p);
