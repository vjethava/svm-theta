function [k]=get_diameter_bound(N, gamma)
if gamma > ((N-2)/(N-1)) 
  k = 1; 
elseif gamma >= 0.5
  k = 2; 
else
  subgraph_size = [2:N]'; 
  inner_fn = @(S) inner(S, gamma); 
  k_bounds = arrayfun(inner_fn,subgraph_size);
  [k, idx] = max(k_bounds); 
end
end

function [k]=inner(N, gamma)
  if gamma > 1/(N-1)
    s = ceil(gamma * (N - 1)) + 1;
    t = floor(N/s);
    switch mod(N, s)
      case 0
        k = 3*t - 3;
      case 1
        k = 3*t - 2;
      otherwise
        k = 3*t - 1;
    end
  else
    k = N - 1; 
  end

end
