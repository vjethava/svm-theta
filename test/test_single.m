%%% This script presents a simple example to check that SVM is working

clc; clear all; close all; 

n = 10; 
p = 0.5; 
warning('on', 'lmb:verbose'); 

%% Graph generation 
A2 = erdosRenyi(n, p); 
[K2, rho2] = getPsdKfromA(A2, 'use_eigen', true); 

%% Using SVM - finds independent set 
model2 = findIndependentSetCSVM(K2); 
alpha2 = model2.x; 

%% Using MKL - finds dense subgraph, so need to complement
G{1} = Graph(A2);
G{1} = G{1}.complement; 
% We add a dummy all ones graph as MKL crashes for wierd reasons!!
G{2} = Graph(ones(n, n)); 
[K, rho, Y, C]=SimpleMKL.get_mkl_kernel(G);
options = SimpleMKL.get_mkl_options();
[alpha, model] = SimpleMKL.run_mkl(K, Y, C, options, 1);

% Results - we should see that the two methods 
% MKL - dense graph over complement 
% SVM - independent set over original using libsvm 
% yield similar results
fprintf(2, '\nSupport vectors using [LIBSVM MKL]\n'); 
disp([alpha alpha2] ) 