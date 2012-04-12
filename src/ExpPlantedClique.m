classdef ExpPlantedClique
% Top-level class for all planted clique experiments (JMLR version)
    properties (Constant)
        logFile = '~/netgem-research/tmp/pcLog_v1.txt'; 
        rawFile = '~/netgem-research/tmp/pcRaw_v1.mat'; 
    end
    methods (Static)
        function k = getK(n, t, p), if nargin < 3, p=0.5; end, k = ceil(2*t*sqrt(n*(1-p)/p)); end
        function rho = getRho(n, k, p), if nargin < 3, p=0.5; end, rho=2*sqrt(n*p*(1-p)) + k*p; end
        function K = getLS(A, rho), n = size(A, 1); K = A/rho + eye(n); end
        function [result, x] = singleRun(n, t, p)
            k = ExpPlantedClique.getK(n, t, p); 
            [~, A] = getPlantedClique(n, 1 - p, k);
            rho = ExpPlantedClique.getRho(n, k, p); 
            K = ExpPlantedClique.getLS(A, rho); 
            [model] = findIndependentSetCSVM(K); 
            result = struct; 
            x = model.x; 
            result.n = n; 
            result.p = p; 
            result.t = t; 
            result.true_k= k; 
            result.v_by_k = model.v / k;
            result.t_bound = 1+1/t; 
            k_stat = PCExperiment.get_k_stat(model, k);
            result = catstruct(result, k_stat); 
        end
        function jmlr_experiment()
            
            warning('on', 'lmb:verbose');
            warning('on', 'lmb:debug'); 
            n = 20000; 
            ns = 1;
            c_vec = [0:5];
            t_vec = [1:0.5:5];
            rfile = ExpPlantedClique.rawFile;
            for i=1:ns
                for j=c_vec
                    for t=t_vec
                        p = 0.5 * n^(-1/3 + j/15); 
                        fprintf(2,'starting single_run() iteration\n'); 
                        [r, x] = ExpPlantedClique.singleRun(n, t, p); 
                        if ~exist(rfile, 'file'), result = r;                             
                        else load(rfile); result(end+1) = r; end
                        save(rfile, 'result'); 
                        clear result;                        
                    end    
                end
            end
        end
        function extract_plot() 
            close all; 
            n = 20000; 
            load(ExpPlantedClique.rawFile); 
            pr = [result.p]';
            tr = [result.t]'; 
            vbkr = [result.v_by_k]'; 
            tbr = [result.t_bound]'; 
            t_vec = [1:0.5:5]; 
            nt = length(t_vec); 
            for c=1:5        
                pc = 0.5 * n^(-1/3  + c/15); 
                vv = zeros(nt, 1); muv = zeros(nt, 1); 
                tt = zeros(nt, 1); mut = zeros(nt, 1); 
                for tc=1:nt
                    t = t_vec(tc);
                    idx = find([tr == t] .* [pr == pc]);
                    vv(tc) = mean(vbkr(idx)); muv(tc) = var(vbkr(idx)); 
                    tt(tc) = mean(tbr(idx)); mut(tc) = var(tbr(idx)); 
                end
                h = figure; 
                plot(t_vec, tt, 'r--'); 
                hold on; box on; 
                errorbar(t_vec, vv, muv, muv, 'b-*'); 
                title(sprintf('p=%g', pc) ); 
                ll = {};
                ll{1} = '1+1/t'; 
                ll{2} = '\omega(K)/\vartheta (G)';                
                legend(ll, 'Location', 'NorthEast'); 
                name = sprintf('~/netgem-research/jmlr/trunk/epc/c-%d' , c);
                saveas(h, name, 'fig'); 
                saveas(gcf, name, 'epsc'); 
            end
            keyboard; 
        end
    end
end
