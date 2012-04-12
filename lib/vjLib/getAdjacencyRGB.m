function [Argb]=getAdjacencyRGB(A, varargin)
%% SHOWADJACENY To use for showing matrix in RGB  format
%
% Optional arguments: 
% -------------------
% 'background': 0-dark (default), 1-light, 2-grayscale
% 'color':      0-green (default), 1-red, 2-blue
%
% Example: 
% --------
% A = full(erdosRenyi(20, 0.3)); 
% Argb = getAdjacencyRGB(A, 'background', 1, 'color', 1);
% imshow(Argb); 
%
[background, varargin]=process_options(varargin, 'background', 0); 
[color, varargin] = process_options(varargin, 'color', 0); 
N = size(A, 1); 
if background==0
    filler = zeros(N, N); 
    data = A; 
elseif background == 1 
    filler = ones(N, N); 
    data = 1 - sign(A); 
else
    filler = A; 
    data = A; 
end

if color==0 % green
    Argb = cat(3, filler, data,  filler);
elseif color==1 % red
    Argb = cat(3, data, filler,  filler); 
else % blue
    Argb = cat(3, filler,  filler, data); 
end
