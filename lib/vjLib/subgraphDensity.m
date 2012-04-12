function [density]=subgraphDensity(A, node_set)
%% SUBGRAPHDENSITY Returns the density for the given node_set in
%% graph with adjacency matrix A (assumed symmetric); 
% 
%  Usage: [density] = subgraphDensity(A, node_set)
%  
    As = A(node_set, node_set); 
    num_edges = nnz(As); 
    num_nodes = size(As, 1); 
    density = num_edges/(num_nodes *(num_nodes - 1)); 
end

