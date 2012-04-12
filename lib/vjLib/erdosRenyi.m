function G = erdosRenyi (n, p)
% ERDOSRENYI Generate the erdos-renyi random graph G(n, p) 
%
% Usage: G = erdosRenyi (n, p)
%
% Returns
% -------
% G: the adjacency matrix for the generated graph
%
% Expects
% -------
% n: number of vertices in the graph 
% p: the edge probability
%
G = rand(n, n) < p; 
G = triu(G, 1);
G = G + G';
% G = sparse(G); 