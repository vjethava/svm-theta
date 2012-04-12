function [X_2, already_found, ret_val] = recursive_mine(obj, X, p_x)
X_2 = [];
already_found = 0; 
ret_val = 0;
% algorithm_variables
augmented = 1; 

while augmented
% step 6-10
[p_x2, is_valid] = obj.get_p_x(X, p_x); 
if ~is_valid
  return 
end
% step 11
p_x = p_x2;
if obj.check_cgqc(p_x) % && (~obj.check_found(p_x))
  X_2 = p_x;
  ret_val = 1;
  already_found = 0;
  obj = obj.add_clique(X_2);
  return
end
% step 12
[X, augmented] = obj.parent_equiv_sub(X, p_x);
end
% step 13
already_found = 0;
% ngbr_vec = obj.get_sorted_by_theta(setdiff(p_x, X));
max_X = max(X);
ngbr_vec = sort(p_x(p_x > max_X), 'ascend'); 
for ngbr = ngbr_vec
  [x_n, af_n, rv_n] = obj.recursive_mine(union(X, ngbr), p_x);
  if rv_n == 1
    already_found = 1;
  end
end
if already_found == 1
  ret_val = 1;
  return 
end
is_cgqc = obj.check_cgqc(X); 
if is_cgqc 
  already_found = obj.check_found(X); 
  if ~already_found
    ret_val = 1;
    obj = obj.add_clique(X); 
    return
  end
end

  