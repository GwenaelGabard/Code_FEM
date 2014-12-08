function [p,inv_p] = reorder_node(element,nn)
%
% Calculate an optimized numbering of the nodes that minimizes the
% bandwidth of the node graph.
% Usage:
%   [p,inv_p] = reorder_element(element,nn)
% element is a matrix where each column corresponds to an element.
% nn is the number of nodes per element.
% p is a permutation vector.
% p_inv is the inverse permutation vector.

element = element(1:nn,:);
ii = zeros(nn^2,size(element,2));
jj = zeros(nn^2,size(element,2));

for n=1:nn
    ii((1:nn)+(n-1)*nn,:) = element;
    jj((1:nn)+(n-1)*nn,:) = ones(nn,1)*element(n,:);
end
% Build a spare matrix with 1 on the component a_ij when the nodes i and j
% are connected.
K = sparse(ii(:),jj(:),ones(size(element,2)*nn^2,1));
% Use the symetric reversed Cuthill-McKee method to reduce the bandwidth
p = symrcm(K);
% Get the inverse permutation
[~,inv_p] = sort(p);
