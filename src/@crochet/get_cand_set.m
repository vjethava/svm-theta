function [cand_set]=get_cand_set(obj, node_set, p_x)
cand_set = p_x; 
for node = node_set
  for i=1:obj.M
    nbd_i = obj.get_nbd(node, i, p_x); 
    cand_set = intersect(cand_set, nbd_i);
  end
end
%%% Not needed due to intersection property. 
% cand_set = setdiff(cand_set, node_set);
end
