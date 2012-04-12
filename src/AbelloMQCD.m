classdef AbelloMQCD
  %AbelloMQCD Quasiclique detection using Abello'02 algorithm 
  properties
    G
    N
    gamma
    alpha 
  end
  methods
    function obj=AbelloMQCD(G)
      obj.G = G;
      obj.N = size(G, 1); 
      obj.alpha = 0.9; 
    end
    function [phi]=getPhiR(obj, R, gamma)
      if nargin < 3
        gamma = obj.gamma;
      end
      Gr = obj.G(R, R); 
      nR = length(R); 
      phi = nnz(Gr) - gamma*nR*(nR-1)/2;  
    end
    function [phi]=getPhiSR(obj, S, R, gamma)
      if nargin < 4
        gamma = obj.gamma;
      end
      SR = [S R]; 
      phi = obj.getPhiR(SR, gamma); 
    end
    function [nbd]=getNbdGamma(obj, S, gamma)
      T = setdiff([1:obj.N], S); 
      gamma_fn = @(x) obj.getPhiR([S x], gamma); 
      gamma_vec = arrayfun(gamma_fn, T, 'UniformOutput', true); 
      nbd = find(gamma_vec >= 0); 
    end
    function [S_1, gamma_1]=grasp_k(obj, K, max_iter)
      gamma_1 = 0;
      S_1 = []; 
      for i=1:max_iter
        S = obj.construct_dsubgK(K);
        gamma = obj.subgraph_density(S); 
        fprintf(2, 'Initial gamma: %g. Starting local_opt\n', gamma); 
        S = obj.localK_11(S); 
        
        gamma = obj.subgraph_density( S); 
        fprintf(2, '\nFinished Local - gamma: %g.\n', gamma); 
        if gamma == 1
          return 
        end
        disp([i gamma gamma_1]); 
        if gamma > gamma_1
          gamma_1 = gamma; 
          S_1 = S; 
        end
      end
    end
    function [delta]=getDeltaSx(obj, S, x, gamma)
      nbd_x = obj.getNbdGamma([x], gamma); 
      nbd_s = obj.getNbdGamma(S, gamma); 
      dxs = nnz(obj.G(x, S)); 
      nx = length(nbd_x); 
      ns = length(nbd_s); 
      delta = nx + ns*(dxs - gamma*(length(S) + 1));
    end
    function [S]=construct_dsubgK(obj, K)
      gamma_1 = 1; 
      deg_x = sum(obj.G, 1); 
      x= obj.get_rcl([1:obj.N], deg_x);

      S_1 = [x]; 
      S = []; 
      deg_fn = @(x) nnz(obj.G(x, S)); 
      % while (gamma_1 >= gamma) && (length(S) < K)
      while  (length(S) < K)
        S = S_1;
        nbd_s = obj.getNbdGamma(S, gamma_1); 
        nbd_s = setdiff(nbd_s, S); 
        if ~isempty(nbd_s)
          delta_fn = @(x) obj.getDeltaSx(S, x, gamma_1);
          g_vec = arrayfun(deg_fn, nbd_s, 'UniformOutput', true); 
          % g_vec = arrayfun(delta_fn, nbd_s, 'UniformOutput', true); 
          % x = randsample(nbd_s, 1, true, p_vec); 
          x = obj.get_rcl(nbd_s, g_vec); 
    
        else
          s_bar = setdiff([1:obj.N], S);
          nbd_idx = sum(obj.G(S, s_bar)) > 0; 
          nbd_1 = s_bar(nbd_idx); 
          if isempty(nbd_1) 
            return ; 
          end
          deg_vec = arrayfun(deg_fn, s_bar, 'UniformOutput', true); 
          x = obj.get_rcl(s_bar, deg_vec); 
        end
        S_1 = [S x]; 
        gamma_1 = obj.subgraph_density( S_1); 
      end
    end
    function [gamma]=subgraph_density(obj, S)
      ns = length(S); 
      gamma = nnz(obj.G(S, S))/(ns*(ns-1));
    end
    function [S_1]=localK_11(obj, S)
      ns = length(S); 
      S_1 = S; 
      found_better = true;
      count = 1; 
      count_1 = 0; 
      while found_better || ((count < ns) && (count_1 < 10*obj.N)) 
      % while found_better || ((count < ns) )
        fprintf(1, '.'); 
        S = S_1; 
        x = S(count); 
        gamma_s = obj.subgraph_density( S);
        n_bar = setdiff([1:obj.N], S); 
        nsx = nnz(obj.G(x, S));
        y_fn = @(y) ((nnz(obj.G(y, S))-obj.G(y, x)) > nsx);
        y_vec = arrayfun(y_fn, n_bar, 'UniformOutput', true);
        yy = n_bar(y_vec); 

        if ~isempty(yy)
          found_better = true; 
          count = 1; 
          S_1 = setdiff(S, [x]) ;
          deg_fn = @(z) nnz(obj.G(z, S)); 
          hh = @(z) obj.subgraph_density([S_1 z]); 
          gh = arrayfun(hh, yy, 'UniformOutput', true); 
          py = arrayfun(deg_fn, yy, 'UniformOutput', true);
          ys = obj.get_rcl(yy, py); 
          S_1 = [S_1 ys];
          gamma_11 = obj.subgraph_density(S_1); 
          fprintf(2, '\nLocal update - prev: %g gamma: %g iter: %d\n', ...
            gamma_s, gamma_11, count_1); 
        else
          found_better = false;
          count = count + 1; 
        end
        count_1 = count_1 + 1; 
      end
      
    end
    
    function [x, rcl]=get_rcl(obj, idx_vec, p_vec)
      p_min = min(p_vec);
      p_max = max(p_vec); 
      threshold = obj.alpha * p_max + (1 - obj.alpha) * p_min; 
      rcl_idx = p_vec >= threshold; 
      rcl = idx_vec(rcl_idx); 
      if length(rcl) == 1
        x = rcl; 
      else
        x = randsample(rcl, 1);
      end
    end
  end
  
end

