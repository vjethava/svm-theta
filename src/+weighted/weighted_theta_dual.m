function [dual_obj, result, t_cvx] = weighted_theta_dual(S, ...
                                                  require_equality, require_acute, show_verbose)
% [H , dual_val, t_cvx] = weighted_theta_dual(S,require_equal)
% 
% Computes dual of the weighted theta as:
%
% max_{X, t}    -t
% s.t.          X_{ii}  =   1; 
%               X_{ij} <=   S_{ij}; 
%               X_{(n+1)i} >-= t; 
% 
% Given by: 
%
% 

tb_cvx = tic; 
n = size(S, 1); 
if nargin < 3 
    require_acute = false; 
end

if nargin < 4 
    show_verbose = true; 
end

bZ = zeros(n, 1); 

cvx_quiet (~show_verbose); 
cvx_solver sdpt3; 
cvx_begin sdp
    variable Y(n, n) symmetric; 
    variable Z(n, n) symmetric; 
    if require_acute
        variable mu0(n, n); 
    end
    variable beta_1(n, 1); 
    variable phi; 
    variable muW(n, n); 
    variable xi(n, 1); 
    variable dual_obj ; 
    maximize dual_obj; 
    subject to
        dual_obj ==  - (sum(beta_1, 1) + phi + trace(muW * S) ); 
        for i=1:n
            for j=1:n   
                if require_acute
                    mu0(i, j) >= 0; 
                end
                muW(i, j) >= 0;                
            end
            if require_acute
                mu0(i, i) == 0; 
            end
            muW(i, i) == 0; 
        end
        if (~require_equality)
            for i=1:n
                xi(i) <= 0;
            end
        end
        
        sum(xi) == -0.5; 
        [Y xi; xi' phi] == semidefinite(n+1); 
        if require_acute
           Z ==  - Y + diag(beta_1) - mu0 + muW ;
        else
           Z ==  - Y + diag(beta_1)  + muW;
        end
        Z == semidefinite(n); 
cvx_end
result=struct('Y', Y, 'dual_obj', dual_obj,  'muW', muW, ...
              'xi', xi, 'phi', phi, 'beta', beta_1); 
H = diag(beta_1) + muW; 
result.H =  H;


w_theta = 1/dual_obj^2; 
if show_verbose
    fprintf(2, 'dual_obj: %g w_theta: %g\n' , dual_obj, w_theta);


if exist('mu0', 'var') 
    fprintf(2, 'max(abs(mu0)): %g\n' , max(max(abs(mu0))));
end
fprintf(2, 'trace(Z): %g phi: %g total: %g\n', trace(H), phi, trace(H)+phi);end  

t_cvx = toc(tb_cvx); 
