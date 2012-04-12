function [model] = findIndependentSetCSVM(K, varargin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% FINDINDEPENDENTSETCSVM  Runs the C-SVM formulation to obtain the 
% independent set solving the problem 
%       \max_{x >= 0} 2e'x-x'Kx
%
% Usage: [model] = findIndependentSetCSVM(K, varargin)                   
%                                                                        
% Expects:                                                               
% ----------------------------------------                               
% K               :                                                      
% Optional arguments
% 'rho'           : Factor for psd 
% 'c_threshold'   : Threshold for selecting quasiclique (Default: c_\rho)
% 's_threshold'   : If specified, returns the quasi-clique formed
%                   by top-s support vectors             
% 'tmp_dir'       : Temporary directory 
% Returns:                                                               
% ----------------------------------------                               
% model           : Model file containing SVM results and post-processing      
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

%% Pre-processsing

% add path for LIBSVM
if exist('process_options') == 0
  run('../startup.m'); 
end

if exist('svmtrain')==0
    warning('lmb:debug', 'LIBSVM path is hard coded'); 
    LIBSVM_ROOT='~/netgem-research/lib/libsvm-3.1'; 
    [svmpath, varargin] = process_options(varargin, 'svmpath', ...
        [LIBSVM_ROOT '/matlab']); 
    addpath(svmpath); 
end
if exist('svmWriteKernel') == 0
    addpath('~/netgem-research/cvx_work/vjLibsvmApi'); 
end
obj = @(x, K) 2*sum(x) - x'*K *x; 

%% Initialization 
assert(size(K, 1) == size(K, 2)); 
N = size(K, 1); 
[C, varargin] = process_options(varargin, 'C', N); 
% avoid cluttering working directory 
% if exist('TMP_DIR', 'var')
%     [tmp_dir, varargin]=process_options(varargin, 'tmp_dir', TMP_DIR);     
% else
%      tmp_dir = '~/netgem-research/tmp';
% end

[kernelFile, varargin] = process_options(varargin, 'kernelFile', tempname);
[modelFile, varargin] = process_options(varargin, 'modelFile', tempname); 
[rho, varargin] = process_options(varargin, 'rho', 0.0); % unspecified 
[c_threshold, varargin] = process_options(varargin, 'c_threshold', 0); %unspecified
[gamma_req, varargin] = process_options(varargin, 'gamma_req', 0); 
[SVMROOT, varargin] = process_options(varargin,  'SVMROOT', '~/netgem-research/lib/libsvm-3.1');
[plot_gamma, varargin] = process_options(varargin, 'plot_gamma', false); 
[gamma_file, varargin] = process_options(varargin, 'gamma_file', tempname); 

%% Writing K to kernelFile
tic;
% augment the k-matrix to put zero as alternate class
K_augment = [2*K zeros(N, 1); zeros(1, N) 0]; 
svmWriteKernel(K_augment, kernelFile, [ones(N, 1) ; -1]); 
tb = toc; 
warning('lmb:verbose', 'Time taken for writing file: %g\n', tb); 

%% Running SVM externally
tic; 
unix(sprintf('%s/svm-train -q -s 0 -c %g -t 4 %s %s', SVMROOT, C, kernelFile, modelFile));
tc = toc;
warning('lmb:verbose', 'Time taken by libsvm solver: %g\n', tc); 
%% SVM code
model = svmReadModel(modelFile);
model.kernel = kernelFile; 
model.file = modelFile; % keep the temporary file name used 
b = zeros(N+1, 1); 
b(model.SVs) = model.sv_coef; 
model.x = b(1:N); % note the N-1 as last element is of opposite sign.  
model.t = tc; 
model.rho = rho; 
model.hardDecision = (find(model.x == 1.0))'; 
nHD = length(model.hardDecision); 


model.v = sum(model.x); % obj(model.x, K); % the objective function 
[s_threshold, varargin] = process_options(varargin, 's_threshold', ...
                                          floor(model.v) );

%% Post-processing
[x_sorted, sort_idx] = sort(model.x , 1, 'descend'); % sorted-list
lidx = [1:N]';
% The following function (f1) computes the RHS in alpha <= 1/(1 +
% t1) where t1 is term depending on gamma, s, rho, and n1 in
% equation 
t1 = @(s) 1.0/(1.0 + (1 - gamma_req)*s*(s-1)/(rho*(s - nHD)));
t2 = @(s) ((s==1) || (s <= nHD)); % check that s is valid
f1 = @(s) max(t1(s), t2(s)); 
% This checks that over the fractional support vector
% coefficients. Nothing can be said about the coefficients
% which are zero.
f2 = @(alpha, s) (alpha > 0) && (alpha <= f1(s)); 
% This extracts the adjacency matrix corresponding to subset of
% nodes given by sd. 
subgraph = @(sd) sign(K(sd, sd) - eye(length(sd))); 
% This computes the gamma_exp for given subset
gammaForSubgraph = @(sd) max((length(sd) == 1), ...
                      1.0 - nnz(subgraph(sd))/(length(sd) * (length(sd) - 1)) ...
                      ); 
% This computes gamma_exp for the top-s support vectors
gammaExpFn = @(s) gammaForSubgraph(sort_idx(1:s)); 
% This computes gamma_pred for top-s support vectors
gammaPredFn = @(s) 1.0 - (s - nHD)*rho*(1-x_sorted(s))/(s*(s-1)*x_sorted(s));
    
if plot_gamma && (rho  > 0) 
    
    warning('lmb:verbose', 'plot_gamma starting'); 
    tic; 
    gammaExp = zeros(N, 1); 
    gammaPred = arrayfun(gammaPredFn, lidx);
    prev_density = 1;
    prev_edges = 0; 
    for i=1:N
          if (i==1) || (i <= nHD) 
              curr_density = NaN; 
              curr_edges = 0; 
          else 
            curr_row = sort_idx(i);
            curr_col = sort_idx(1:(i-1));
            curr_edges = prev_edges + 2*nnz(K(curr_row, curr_col)); 
            curr_density = 1.0 - curr_edges / (i *(i-1));
          end
        prev_edges = curr_edges; 
        prev_density = curr_density; 
        gammaExp(i) = curr_density; 
    end
    gammaPred(isinf(gammaPred)) = NaN; 
    gammaPred(gammaPred<0) = NaN; 
    tg = toc;     
    warning('lmb:verbose', 'plot_gamma finished in %s. Plotting', tg); 
    h = figure; 
    semilogx(lidx, gammaExp, vjGetLineStyle(0) ); 
    hold on; 
    semilogx(lidx, gammaPred, vjGetLineStyle(1) ); 
    % xlabel('Quasi-clique size'); 
    grid on;
    set(gca, 'FontSize', 20); 
    % legend({'\gamma^{exp}'; '\gamma^{pred}'}, 'Location', 'SouthEast'); 
    plot_exp = [gamma_file '_gamma']; 
    saveas(h, plot_exp, 'epsc'); 
    close(h); 
end

if c_threshold > 0.0 % soft-thresholding 
    sd = find(model.x >= c_threshold); 
    nSD = length(sd); 
    model.softDecision = sd; 
    model.c_threshold = c_threshold; 
    model.gamma_exp = gammaForSubgraph(sd); 
    if rho > 0
        %%% Note: Generic weighted edges - computation of rho
        % Wsd = K(sd, sd) - eye(nSD);
        % [Isd, Jsd, Vsd] = find(Wsd);
        % model.w_min = min(Vsd);
        % rho = 1.0/model.w_min; 
        g_pred = 1 - (nSD - nHD)/(nSD*nSD - nSD)*rho*(1-c)/c;
        model.gamma_pred = g_pred;       
    end
elseif (gamma_req > 0) 
    max_n = floor(log(N)*10); 
    curr_n = 1; 
    above_req = true; 
    max_idx = 2; 
    max_gamma = 0.0; 
    curr_gamma = 0.0; 
    while ((above_req) || (curr_n < max_n)) && (curr_n < N)
        prev_gamma = curr_gamma; 
        curr_n = curr_n + 1; 
        curr_sd = sort_idx(1:curr_n); 
        curr_density = subgraphDensity(K, curr_sd); 
        curr_gamma = 1.0 - curr_density; 
        above_req = (gamma_req <= curr_gamma); 
        mesg =sprintf('n: %d gammaReq: %g gammaExp: %g maxN: %d', curr_n,...
                      gamma_req, curr_gamma, max_n); 
        warning('lmb:debug', mesg); 
        if (curr_gamma > max_gamma) || (curr_gamma > gamma_req)
            max_idx = curr_n;
            max_gamma = curr_gamma; 
        end

    end
    if prev_gamma < gamma_req 
        model.gamma_exp = max_gamma;
        model.softDecision = sort_idx(1:max_idx); 
    else
        model.softDecision = curr_sd; 
        model.gamma_exp = curr_gamma;
    end
    model.gamma_req = gamma_req; 
    
elseif (s_threshold > 0) % compute the top-s support vector result
                         % and report the results.  
    sd = sort_idx(1:s_threshold);  nSD = length(sd);  
    model.gamma_exp = 1 - ...
        (nnz(K(sd, sd) - eye(nSD)))/(nSD * (nSD - 1));
    if (rho > 0)
        model.gamma_pred = gammaPredFn(s_threshold); 
    end
    model.softDecision = sd; 
    model.alpha_s = x_sorted(s_threshold); 
    model.s_threshold = s_threshold; 
end
%% Clean-up
% unix(sprintf('rm %s %s', kernelFile, modelFile));
delete(kernelFile); 
delete(modelFile); 
end

function [y] = neglog(x)
    if x >= 0
        y = x;
    else
        y = log(-x);
    end
end

