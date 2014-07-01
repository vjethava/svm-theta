n = 10000; 
p = 0.5; 
A = getPlantedClique(n, p, ceil(6*sqrt(n))); 
K = getPsdKfromA(A); 
f = -ones(n, 1); 
lb = zeros(n, 1);
ub = ones(n, 1);
tic; [model] = findIndependentSetCSVM(K); t_svm = toc 
x_svm = model.x;

% tic; 
% [x_qp, f_qp] = quadprog(K, f, [], [], [], [], lb, ub); 
% t_qp = toc 

tic; 
[x_admm, f_admm] = quadprog_admm(K, f, 0 , lb, ub, 1.0, 1.0); 
t_admm = toc 

[norm(x_svm - x_admm,1) ]
