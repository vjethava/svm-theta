function [S] = ananthram_generator(n, k, epsilon)
% 
%   Generates weighted graph with clique of size k 
t = floor(n/k); 
A =  kron( eye(t), ones(k, k));
S = (epsilon)*(ones(n, n) - A) + (1 - epsilon)*A; 
for i=1:n
    S(i, i) = 0; 
end
