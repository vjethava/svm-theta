classdef ExpMultDimacs 
% Experiment which compares density obtained for union graph vs MKL-1 
% implementation for Support Vectors 
    properties (Constant)
        res_dir = '~/netgem-research/jmlr/trunk/mult-dimacs';
        res_file = [ExpMultDimacs.res_dir '/result-v1.mat'];
        fig_dir = ['~/netgem-research/jmlr/trunk/t-fig']; 
    end
    
    methods 
        function obj = ExpMultDimacs()
            set(0,'defaultaxesfontsize',15);
            set(0,'defaulttextfontsize',15);
            % We want the line widths to be 2 points instead of the usual 1.
            set(0,'defaultaxeslinewidth',1);
            set(0,'defaultlinelinewidth',1);
            set(0, 'DefaultLineMarkerSize', 4, ...
                    'DefaultLineLineWidth', 2);
            if exist(obj.res_dir, 'dir') == 0
                mkdir(obj.res_dir);
            end 
            result = zeros(14, 6);
            alpha = cell(14, 1);  
            G = cell(14, 1);
            rho = cell(14, 1); 
            for i=1:14
                [r, l, alpha_c, G_c, rho_c] = ExpMultDimacs.run_graph(i); 
                result(i, :) = ExpMultDimacs.get_table_row(l); 
                alpha{i} = alpha_c;
                rho{i} = rho_c; 
                G{i} = G_c; 
                ExpMultDimacs.plot_Tc(r, G_c, i, alpha_c); 
                % last_stats{i} = l; 
            end
            b = [result(:, 2), result(:, 3) result(:, [4 5 6])];
            b_file = [ExpMultDimacs.res_dir '/table1.tex'];
            matrix2latex(b, b_file); 
            save(ExpMultDimacs.res_file, 'result', 'alpha', 'G', 'rho');     
        end
        
    end
    methods (Static) 
        function plot_Tc(r, G, i, alpha)
            gamma = {r.mkl}; 
            names = {r.graph};
            m = length(gamma); 
            nt = size(gamma{1}, 1); 
            f = figure; 
            box on; grid off ; hold on; 
            
            legends  = cell(m, 1); 
            % legends{1} = '\gamma_r = 1'; 
            
            for l=1:m
                c_gamma = gamma{l}(:, 2)/G{l}.density; 
                plot([1:nt], c_gamma, vjGetLineStyle(l+1)); 
                c_name =  G{l}.name;
                c_name = regexprep(c_name, '\_complement', '^c'); 
                c_name = regexprep(c_name, '\_', '\\\_'); 
                legends{l} = c_name; 
            end
            legend(legends, 'Location', 'NorthEast');
            nSV = nnz(alpha);
            plot([1:nSV], ones(1, nSV), 'k--');
            plot(0, 0, 'k.'); 
            fig_name = [ExpMultDimacs.fig_dir sprintf('/exp-%d', i)];
            saveas(f, fig_name, 'fig'); 
            % saveas(f, fig_name, 'epsc'); 
        end
        function [gammaRow] = get_table_row(r)
            densityRatio =  r.density(2, :)./r.density(3, :); 
            minDensity = min(densityRatio); 
            maxDensity = max(densityRatio); 
            avgDensity = mean(densityRatio); 
            gammaRow = [r.nSV(1) r.nSV(2) r.nSV(3) maxDensity avgDensity minDensity]; 
        end
        function [result, finalStat, mkl_alpha, G, rho]=run_graph(graph_family_idx)              
            sparse_check = true; 
            % dimacs graphs
            names = DIMACS.graph_names{graph_family_idx}; 
            G = DIMACS.get_graphs(graph_family_idx, sparse_check);
            mkl_opts = SimpleMKL.get_mkl_options(); 
            [K, rho, Y, C] = SimpleMKL.get_mkl_kernel(G);
            [mkl_alpha, mkl_model] = SimpleMKL.run_mkl(K, Y, C, mkl_opts);  
            [mkl_T, mkl_alpha_s] = SimpleMKL.postprocess(K, mkl_alpha, mkl_model, rho); 
            Gu = Graph.union_graph(G); 
            Ku = getPsdKfromA(Gu, 'use_eigen', 1); 
            [union_model] = findIndependentSetCSVM(Ku); 
            % collecting stats. 
            f = figure;
            m = length(G); 
            mkl_res =cell(m, 1); 
            union_res = cell(m, 1); 
            final_density =  zeros(3, m);
            for i=1:m
                final_density(3, i) = G{i}.density; 
            end
            final_nsv = zeros(1, 3); 
            % number of support vectors for union-graph 
            final_nsv(1)  = union_model.nSV(1); 
            % number of support vectors for mkl-1
            final_nsv(2) = nnz(mkl_alpha_s); 
            final_nsv(3) = nnz(mkl_alpha);
            for i=1:length(G)
                
                subplot(m, 1, i); 
                hold on; 
                grid on; 
                set(gca, 'FontSize', 10); 
                title(strrep(G{i}.name, '_', '\_')); 
                [mo, mg, ma] = Graph.get_ordered_density(G{i}, mkl_alpha_s); 
                [uo, ug, ua] = Graph.get_ordered_density(G{i}, union_model.x); 
                
                msv = [1:length(mg)]; 
                semilogx(msv, mg, 'b-*', 'LineWidth', 2, 'MarkerSize', 3); 
                final_density(2, i) = mg(end); 
                usv = [1:length(ug)]; 
                semilogx(usv, ug, 'r--d', 'LineWidth', 2, 'MarkerSize', 3); 
                final_density(1, i) = ug(end); 
                mkl_res{i} = [mo mg ma]; 
                union_res{i} = [uo ug ua]; 
            end
            result = struct('graph', names , 'mkl', mkl_res, ...
                'union', union_res); 
            fileName = sprintf('%s/%d.eps', ExpMultDimacs.res_dir, ...
                graph_family_idx);
            finalStat =  struct('density', final_density, 'nSV', final_nsv); 
            % print(f, '-depsc', fileName);  
            
            close(f);
        end
        
        
        
    end
    
end
