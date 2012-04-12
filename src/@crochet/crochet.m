classdef crochet < handle
  %CROCHET Implements the CROCHET algorithm for finding cross graph
  % quasicliques given in [Pei'2005] 
  properties (SetAccess='private')
    graph_cell 
    gamma
    min_size
    min_sup
  end
  properties 
    maximal_cliques
    nbd_cell
    cand_G
    cand_V
  end
  properties (Dependent)
    M 
    N 
  end
  methods
    function obj=add_clique(obj, X)
      fprintf(2, 'Adding maximal clique\n'); 
      disp(X); 
      obj.maximal_cliques{end+1} = X;
    end
    function [already_found, par_clique]=check_found(obj, node_set)
      check_fn = @(par_set) all(ismember(node_set, par_set));
      found_vec = cellfun(check_fn, obj.maximal_cliques, ...
        'UniformOutput', true);
      already_found = any(found_vec); 
      par_clique = [];
      if already_found
        par_clique = obj.maximal_cliques{found_vec == 1};
        fprintf(2, 'check_found %s in %s\n', ...
          mat2str(node_set) , mat2str(par_clique));
        % disp(par_clique);
      end
    end
    function [m]=get.M(obj), m = length(obj.graph_cell); end
    function [n]=get.N(obj), n = size(obj.graph_cell{1}, 1); end 
    function obj = crochet(ip_graph_cell, ip_gamma, ip_min_size, ip_min_sup) 
      obj.graph_cell = ip_graph_cell; 
      obj.gamma = ip_gamma;
      obj.min_size = ip_min_size; 
      obj.min_sup = ip_min_sup; 
      obj.maximal_cliques = {}; 
      obj.nbd_cell = cell(obj.N, obj.M); 
      [obj.nbd_cell{:}] = deal({}); 
    end
    [clique_cell] = find_maximal_cliques(obj);
    function [nbd_i]=get_nbd(obj, node, i, p_x)
      n_i = length(p_x);
      k_i = obj.get_diameter_bound(n_i, obj.gamma(i));
      [nbd_i] = obj.get_k_nbd(node, k_i, i, p_x);
    end
    [valid_vertex] = check_vertex(obj, node, p_x);
    [cand_set] = get_cand_set(obj, node_set, p_x);
    function [nbd] = get_k_nbd(obj, node, k, i, p_x)
      nbd = [];
      if length(obj.nbd_cell{node, i}) >= k
        nbd = intersect(obj.nbd_cell{node, i}{k}, p_x);
        return
      end 
      old_p_x = p_x; 
      if k==1
        
        p_x = [1:obj.N]; 
        n_idx = obj.graph_cell{i}(node, p_x);
        nbd = p_x(n_idx==1);
      else
        [first_nbd] = obj.get_k_nbd(node, 1, i, p_x);
        for l=first_nbd
          n_l = obj.get_k_nbd(l, k-1, i, p_x);
          n_l = setdiff(n_l, [node]);
          nbd = union(n_l , nbd);
        end
        nbd = union(nbd, first_nbd);
      end
      obj.nbd_cell{node, i}{k} = nbd; 
      nbd = intersect(old_p_x, nbd); 
    end
    function [p_x_sorted]=get_sorted_by_theta(obj, p_x)
      nbd_i = @(i) [p_x(1:(i-1)) p_x((i+1):length(p_x))]; 
      theta_fn = @(i) obj.get_min_theta(p_x(i), nbd_i(i)); 
      theta_vec = arrayfun(theta_fn, [1:length(p_x)]); 
      [theta_sorted, sort_idx] = sort(theta_vec, 'descend'); 
      p_x_sorted = p_x(sort_idx); 
    end
        
    function [min_theta]=get_min_theta(obj, node, p_x)
      theta_fn = @(i) obj.get_theta_i(node, i, p_x); 
      theta_vec = arrayfun(theta_fn, [1:obj.M]); 
      min_theta = min(theta_vec); 
    end
    function [is_cgqc]=check_cgqc(obj, X)
      is_cgqc = true; 
      nX = length(X); 
      if nX < obj.min_size
        is_cgqc = false; 
      else 
        p_i = @(i) [X(1:(i-1)) X((i+1):nX)];
        j = 1;
        while (is_cgqc) && (j <= obj.M)
          deg_fn = @(i) obj.get_degree_i(X(i), j, p_i(i));
          deg_graph_j = arrayfun(deg_fn, [1:nX]) ;
          is_cgqc = all(deg_graph_j  >= obj.gamma(j)*(nX - 1));
          j = j+1;
        end
      end
    end
    function [X_2, augmented]=parent_equiv_sub(obj, X, p_x)
      X_2 = X; 
      nX = length(X);
      n_x = setdiff(p_x, X); 
      augmented = false;
      for u = X
        deg_fn = @(i) obj.get_degree_i(u, i, p_x);
        req_vec = obj.gamma.*(max(obj.min_size, nX)-1); 
        deg_vec = arrayfun(deg_fn, [1:obj.M]); 
        is_in_cgqc = (req_vec == deg_vec); 
        if(any(is_in_cgqc))
          
          g_i_cond = find(is_in_cgqc) ;
          for i = g_i_cond
            nbd_u = obj.get_k_nbd(u, 1, i, n_x);
            X_2 = union(X_2, nbd_u);
          end
        end
      end
      if length(X_2) > length(X)
        fprintf(2, 'Found augment nodes: %s for %s\n', ...
          mat2str(X_2), mat2str(X)); 
        augmented = true;
      end
    end
    [X_2, already_found, ret_val] = recursive_mine(obj, X, p_x);
    function [p_x2, is_valid]=get_p_x(obj, node_set, p_x)
      y = obj.get_cand_set(node_set, p_x); 
      p_x2 = union(node_set, y); 
      val_idx = false;
      % at least one node got thrown away
      while ~all(val_idx)
        prune_fn = @(node) obj.check_vertex(node, p_x2);
        val_idx = arrayfun(prune_fn, p_x2);
        p_x2 = p_x2(val_idx);
      end
      x_in_p = all(ismember(node_set, p_x2));
      is_valid = true; 
      % prune y - and check validity
      if length(p_x2) < obj.min_size
        is_valid = false;
      elseif ~x_in_p
        is_valid = false; 
      else
        already_found = obj.check_found(p_x2); 
        if already_found
          is_valid = false;
        end
      end
    end
  end
  methods (Access='private') 
    [valid_vertex] = check_vertex_i(obj, node, i, p_x);
    function [deg_i]=get_degree_i(obj, node, i, p_x)
      deg_i = nnz(obj.graph_cell{i}(node, p_x)); 
    end
    function [theta] = get_theta_i(obj, node, i, p_x)
      theta = obj.get_degree_i(node, i, p_x)/obj.gamma(i); 
    end
  end
  methods (Static)
    [diameter_bound] = get_diameter_bound(graph_size, gamma);     
  end
end

