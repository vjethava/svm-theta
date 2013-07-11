function [w_theta, dual_obj, xi] = dual_theta_2(S)
%% [dual_obj, xi] = dual_theta_2(S)
% 
% Implements simpler dual theta 
% 
% - \min \sqrt( \xi^\top (I + S) \xi ) 
%   s.t. 
%        \sum_i \xi_i = 1, \xi_i \ge 0 \forall i 
%
% S_ij = similarity between nodes i and j 
% S_ii = 0
% 
n = size(S, 1); 
cvx_quiet true; 

cvx_begin  
    variable xi(n, 1);
    variable H(n, n) symmetric; 
    minimize quad_form( xi, (eye(n) + S)  );
    subject to 
        xi == simplex(n); 
cvx_end

dual_obj = - sqrt( cvx_optval ); 
w_theta = 1 / cvx_optval; 

 
