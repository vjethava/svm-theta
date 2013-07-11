function [w_theta, X , primal_obj , t_cvx] = weighted_theta(S, ...
                                                  require_equal,require_acute, show_verbose)
% [w_theta, X , t, t_cvx] = weighted_theta(S,require_equal)
% 
% Computes weighted theta as:
%
% min_{X, t}    - t
% s.t.          X_{ii}  =   1; 
%               X_{ij} <=   S_{ij}; 
%               X_{(n+1)i} >-= t; 
% 
% Theta is defined as:  w_theta = (1/t)^2 
% 
tb_cvx = tic;    
n = size(S, 1);
tb_cvx = tic;    
if nargin < 4
    show_verbose = true; 
end

if nargin < 3
    require_acute = false; 
end

cvx_quiet (~show_verbose); 
cvx_solver sdpt3; 
cvx_begin sdp 
    variable X(n+1, n+1) symmetric;
    variable t;

    minimize  -t ; 
    subject to 
        X == semidefinite(n+1);         
        for i=1:(n+1)
            X(i, i) == 1;
        end
        for i=1:n-1
            for j=i+1:n
                X(i,j) <=  S(i,j);
                % acute labelling required. 
                if require_acute 
                    X(i, j) >= 0; 
                end
                
            end
        end
        for i=1:n
            if require_equal
                X(n+1,i) == t;
            else
                X(n+1,i) >= t;
            end
        end
cvx_end
primal_obj = cvx_optval; 
w_theta = 1/t^2; 
if show_verbose
    fprintf(2, 'primal_obj: %g w_theta: %g\n' , primal_obj, w_theta);
end
t_cvx = toc(tb_cvx);
