classdef Graph < double
  %GRAPH Generic weighted undirected graph class
  properties (Constant, Access='private')
    DIMACS_DIR='../data/';
    TMP_DIR='/tmp';
  end
  %% Static methods - graph generation 
  methods (Static)
    function [svmStruct]=csvm(K)
        % K = obj.K; 
        N = size(K, 1);
        K_augment = [2*K zeros(N, 1); zeros(1, N) 0]; 
        xx = [1:(N+1)]'; 
        yy = [ones(N, 1); -1];
        kfun = @(xi, xj) K_augment(xi, xj);
        svmStruct = svmtrain(xx, yy, 'Kernel_Function', kfun, 'method', 'SMO'); 
    end
    function [obj]=ErdosRenyi(N, p)
      G = erdosRenyi(N, p); 
      obj=Graph(G); 
      obj.type = 'ErdosRenyi'; 
      obj.name = sprintf('G(%d,%.4g)', N, p); 
      obj.p = p; 
    end
    
    function [node_order, subgraph_density, alpha_sorted] = get_ordered_density(G, alpha, c)
        if nargin < 3
            c = 0; 
        end
        A = double(G); 
        [alpha_sorted, node_order] = sort(alpha, 'descend');
        valid_idx = alpha_sorted > c; 
        node_order = node_order(valid_idx); 
        alpha_sorted = alpha_sorted(valid_idx); 
        ns = length(node_order); 
        subgraph_density = ones(ns, 1);   
        curr_edges = 0; 
        for i=2:ns
            curr_node = node_order(i);
            prev_nodes = node_order(1:(i-1)); 
            edges_added = sum(A(curr_node, prev_nodes)); 
            curr_edges = curr_edges + edges_added;  
            subgraph_density(i) =  2*curr_edges/(i*(i-1)); 
        end
    end
    function [Gu]=union_graph(G)
        num_graphs = length(G); 
        A = double(G{1}); 
        N = G{1}.N
        for j=2:num_graphs
            assert(G{j}.N == N, 'Graphs should be of same size'); 
            A = A + double(G{j});
        end
        Gu = Graph(sign(A)); 
        Gu.name = 'union_graph';
    end
    
    function [obj]=DIMACS(name)
      edge_file = [Graph.DIMACS_DIR '/ascii/' name '.clq'];
      edges = load(edge_file);
      E = size(edges, 1);
      N = max(max(edges));
      A = zeros(N, N, 'double'); 
      for i=1:E   % this part is slow - see bin2asc.c and genbin.h
          A(edges(i, 1), edges(i, 2)) = 1;
      end
      A = A + A'; % symmetric edges ij \in E => ji \in E
      obj = Graph(full(A)); 
      obj.type = 'DIMACS';
      obj.name = name;      
    end
    
    function [obj]=PlantedClique(N, p, k)
      A = getPlantedClique(N, p, k); 
      obj=Graph(A); 
      obj.type = 'PlantedClique';
      obj.p = p;
      obj.k = k; 
      obj.name = sprintf('G(%d,%.4g,%d)', N, p, k); 
    end
    
  end
  %% Properties dependent on the graph 
  properties (Dependent)
      N                % number of nodes
      E                % number of edges
      max_deg          % maximum degree
      avg_deg          % average degree
      density          % graph density 2*E/N*(N-1)
  end
  
  properties (SetAccess='protected')
    type
    p                 % valid only for random graphs
    k                 % clique size
  end
  
  %% SVM properties
  properties
    name
    K 
    rho 
    model
  end
  %% Generic functions
  methods 
    
    function [k, r] = compute_psd(obj)    
        v = eig(double(obj)); 
        r = -v(1);
        k = eye(obj.N, obj.N) + double(obj)/r; 
    end
    function [k]=get.K(obj)
      k = eye(obj.N, obj.N) + double(obj)/obj.rho; 
    end
    function disp_rho(obj)
        fprintf(2, 'rho: %g\n', obj.rho); 
    end
    function [G]=complement(obj)
      n = obj.N; 
      a_c = ones(n, n, 'double') - eye(n, n, 'double') - double(obj); 
      G = Graph(a_c); 
      if ~isempty(obj.name)
          G.name = [obj.name '_complement'];
      end
    end
    function obj=set.K(obj, k)
      obj.K = k;
    end    
    function obj=set.rho(obj, r)
      obj.rho = r; 
    end
    function disp(obj)
      if ~isempty(obj.name)
        fprintf(1, '\tname:\t%s\n', obj.name);
      end
      fprintf(1, '\tG.N:\t%d\n', obj.N); 
      fprintf(1, '\tG.E:\t%d\n', obj.E); 
    end
  end
  
  %% Dependent property methods
  methods
    function obj=Graph(data)
      obj=obj@double(data); 
    end   
    function [n]=get.N(obj)
      n = size(obj, 1);
    end
    function [e]=get.E(obj)
      e = nnz(triu(obj)); 
    end
    function [density] = get.density(obj)
      density = 2*obj.E/(obj.N*(obj.N-1));       
    end
 
%     function [B]=subsref(obj, S)
%       fprintf(2, '%s\n', S.type); 
%       keyboard; 
%       disp(S{:}.subs); 
%     end
    
    function [Delta]=get.max_deg(obj)
      Delta = max(sum(obj)); 
    end
    function [delta]=get.avg_deg(obj)
      if obj.N > 0
        delta = (2*obj.E)/(obj.N);
      else 
        delta = 0; 
      end
    end
  end
  
end

