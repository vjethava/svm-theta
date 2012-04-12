function [valid_vertex]=check_vertex(obj, node, p_x)
valid_vertex= true;
i = 1;
while ((valid_vertex) && (i <= obj.M))
  valid_vertex = obj.check_vertex_i(node, i, p_x);
  i= i+1; 
end
end

