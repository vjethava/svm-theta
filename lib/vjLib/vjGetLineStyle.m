function [lineStyle] = vjGetLineStyle(k)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% VJGETLINESTYLE                                                         
%                                                                        
% Usage: [lineStyle] = vjGetLineStyle(k)                                 
%                                                                        
% Expects:                                                               
% ----------------------------------------                               
% k               : Index                                                      
%                                                                        
% Returns:                                                               
% ----------------------------------------                               
% lineStyle       : Line style string                                                      
%                                                                        
% Note:                                                                  
% ----------------------------------------                               
%                                                                        
% See Also:                                                              
% ----------------------------------------                               
%                                                                        
%                                                                        
% Copyright (c) 2011, Vinay Jethava (vjethava@gmail.com)                 
%                                                                        
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if nargin <  1
    k = 0; 
end
colors = {'b', 'r', 'g' , 'm' , 'c' };
lines = {'-', ':', '-.', '--'};
dots = {'.', 'o', 'x', '+', '*', 's', 'd', 'v', '<', '>', 'p', 'h'};
a = rem(k, length(colors)) + 1;
b = rem(k, length(lines)) + 1;
c = rem(k, length(dots)) + 1;
lineStyle = [colors{a} dots{c} lines{b} ]; 
end
