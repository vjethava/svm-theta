classdef SimpleMKL
    % Refactored class from MKLExperiment2
    methods (Static)
        function [K, rho, Y, C]=get_mkl_kernel(G)
            D = length(G);
            N = G{1}.N;
            K = zeros(N+1, N+1, D);
            rho = zeros(D, 1);
            Y = [ones(N, 1); -1];
            for i=1:D
                A = double(G{i}.complement);
                [ck, crho] = getPsdKfromA(A, 'use_eigen', 1);
                K(:,:, i) = [2*ck zeros(N, 1); zeros(1, N+1)];
                rho(i) = crho;
            end
            C = N;
        end
        function [alpha, model] = run_mkl(K, Y, C, mkl_options, verbose)
            if nargin < 5
                verbose = 1;
            end
            N = size(K, 1);
            D = size(K, 3);
            [beta, w, b, posw, story, objective]=mklsvm(K, Y, C, mkl_options, verbose);
            alpha = zeros(N, 1);
            model = struct('beta', beta, 'w', w, 'posw', posw, ...
                'b' , b, 'story', story, 'obj', objective);
            alpha(posw) = w;
            alpha = alpha(1:(N-1));
        end
        
        function [T_c, alpha_c] = postprocess(K, alpha, model, rho)
            
            M = size(K, 3); 
            T = cell(M, 1); 
            n1 = size(K, 1); 
            n = n1 - 1;
            T_c = [1:n]'; 
            S_c = find(alpha > 0); 
            for l=1:M
                cA = sign(K(:, :, l)); 
                cA = cA(1:n, 1:n); 
                crho = rho(l); 
                T{l} = SimpleMKL.getTfromS(cA, crho, alpha, S_c); 
                T_c = intersect(T{l}, T_c); 
                fprintf(2, '|Sc: %d T{%d}: %d T: %d\n', ...
                    length(S_c), l, length(T{l}), length(T_c)); 
            end
            t_nc = setdiff([1:n]', T_c); 
            alpha_c = alpha;
            alpha_c(t_nc) = 0.0;

        end
        
        
        function [T_l] = getTfromS(A, rho, alpha, S_c)
            if nargin < 4 
                S_c = find(alpha > 0.0); 
            end
            N = size(A, 1); 
            alpha_min = min(alpha(S_c)); 
            %  d_i = sum(A(:, S_c), 2);
            alpha_2 = A(S_c, S_c) * alpha(S_c);
            val = (1 - alpha_min - alpha_2/rho); 
            T_l = S_c(find(val > 0)); 
        end
        
        function [mkl_options]=get_mkl_options()
            mkl_options = struct;
            mkl_options.algo='svmclass'; % Choice of algorithm in mklsvm can be either
            % 'svmclass' or 'svmreg'
            %------------------------------------------------------
            % choosing the stopping criterion
            %------------------------------------------------------
            mkl_options.stopvariation=0; % use variation of weights for stopping criterion
            mkl_options.stopKKT= 0;       % set to 1 if you use KKTcondition for stopping criterion
            mkl_options.stopdualitygap=1; % set to 1 for using duality gap for stopping criterion
            
            %------------------------------------------------------
            % choosing the stopping criterion value
            %------------------------------------------------------
            mkl_options.seuildiffsigma=1e-2;        % stopping criterion for weight variation
            mkl_options.seuildiffconstraint=0.1;    % stopping criterion for KKT
            mkl_options.seuildualitygap=0.01;       % stopping criterion for duality gap
            
            %------------------------------------------------------
            % Setting some numerical parameters
            %------------------------------------------------------
            mkl_options.goldensearch_deltmax=1e-1; % initial precision of golden section search
            mkl_options.numericalprecision=1e-8;   % numerical precision weights below this value
            % are set to zero
            mkl_options.lambdareg = 1e-8;          % ridge added to kernel matrix
            
            %------------------------------------------------------
            % some algorithms paramaters
            %------------------------------------------------------
            mkl_options.firstbasevariable='fullrandom'; % tie breaking method for choosing the base
            % variable in the reduced gradient method
            mkl_options.nbitermax=500;             % maximal number of iteration
            mkl_options.seuil=0;                   % forcing to zero weights lower than this
            mkl_options.seuilitermax=10;           % value, for iterations lower than this one
            
            mkl_options.miniter=0;                 % minimal number of iterations
            mkl_options.verbosesvm=0;              % verbosity of inner svm algorithm
            mkl_options.efficientkernel=0;         % use efficient storage of kernels
        end
    end
end
