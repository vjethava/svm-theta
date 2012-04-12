function [valid_vertex]=check_vertex_i(obj, node, i, p_x)
valid_vertex = true; 
deg_i = obj.get_degree_i(node, i, p_x);
if deg_i < obj.gamma(i) * (obj.min_size - 1)
  valid_vertex=false;
else
  nbd_i = obj.get_nbd(node, i, p_x);
  if length(nbd_i) < obj.min_size - 1
    valid_vertex = false;
  end
end
end