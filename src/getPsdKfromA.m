function [K, rho, ta] = getPsdKfromA(A, varargin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% GETPSDKFROMA Compute K = I + A/(-\lambda_\min(A)) s.t. K is positive 
% semi-definite                                                           
%                                                                        
% Usage: [K] = getPsdKfromA(A)                                           
%                                                                        
% Expects:                                                               
% ----------------------------------------                               
% A               :  Matrix                                                    
%                                                                        
% Returns:                                                               
% ----------------------------------------                               
% K               :  Positive semi-definite matrix                                                    
%                                                                        
% Note:                                                                  
% ----------------------------------------                               
%                                                                        
% See Also:                                                              
% ----------------------------------------                               
%                                                                        
% Copyright (c) 2011, Vinay Jethava (vjethava@gmail.com)                 
%                                                                        
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
[maxN, varargin] = process_options(varargin, 'maxN', 1000); 
[use_max_deg, varargin] = process_options(varargin, 'use_max_deg', ...
                                          0);  
[use_eigen, varargin] = process_options(varargin, 'use_eigen', 0);  
assert(use_eigen * use_max_deg == 0, 'Cannot do both max_deg and eigen'); 
getMaxDeg = @(A) max(sum(sign(A), 1)); 
N = size(A, 1); N = size(A, 1); 
tic; 
if ((N > maxN)  || (use_max_deg == 1)) && (use_eigen == 0)
  rho = getMaxDeg(A);
  method='max_deg';
else
  method='eigen';
  [U, D] = eig(A);
  rho = -min(diag(D));        % lambdaMinOfA < 0 as Tr(A)=0
end
K  = eye(N) + A/rho; % positive semi-definite 
ta = toc; 
warning('lmb:verbose', 'getPsdKfromA() Used: %s Time taken: %g\n', method, ta); 