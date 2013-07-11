% This function checks the dual formulation for weighted theta. 
clear all; close all; clc; 
%%% parameters
n = 20; 
k = 4; 
epsilon = 0.001; 
require_equality = false;
require_acute = true;  
%%% Code generation 
S = weighted.ananthram_generator(n, k , epsilon); 
p = 0.0; 
for i=1:n
    for j=(i+1):n
        if rand() < p
            S(i, j) = rand(); 
            S(j, i) = S(i, j); 
        end
    end
end

show_verbose = false; 

require_equality = true; 
require_acute = true; 

for i=1:4

    require_equality = rem(i, 2); 
    require_acute = floor( i / 2); 
    
[t2, obj2, xi] = weighted.dual_theta_2(S); 
[w_theta, X, primal_obj,  t_cvx] = weighted.weighted_theta(S, require_equality, require_acute, show_verbose);
[dual_obj, res_dual] = weighted.weighted_theta_dual(S, require_equality, require_acute, show_verbose);
fprintf(2, ['req_equal: %d req_acute: %d primal: %g dual: %g w_theta: ' ...
'%g dual_2: %g\n'], require_equality, require_acute, primal_obj, ...
        dual_obj, w_theta, t2); 

end