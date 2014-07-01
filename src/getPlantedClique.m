function [G, G_bar]=getPlantedClique(N, p, k)
%% Generating a planted clique
% 
% [G, G_bar, K, rho]=getPlantedClique(N, p, k)
% 
% N - size of the graph
% p - edge probability 
% k - planted clique size
%
% G - graph with planted clique at first k positions
% G_bar - complement of G
%
if nargin <  1
    varargin = {}; 
    % plant a clique of size 10 in the first 10 positions
    [N, varargin] = process_options(varargin, 'N', 1000); 
    [p, varargin] = process_options(varargin, 'p', 0.5); 
    [k, varargin] = process_options(varargin, 'k', 30); 
end
% G = erdosRenyi(N, p);          
G = rand(N, N) < p; 
G(1:k, 1:k) = 1;                    % planted clique of size k
G = triu(G, 1); 
G = G + G'; 
G_bar = 1-eye(N) - G;                   % complement graph 
% [K, rho] = getPsdKfromA(G, 'maxN', 1000); 
