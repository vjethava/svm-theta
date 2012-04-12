function [prevStream] = vjSetRNG(varargin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% VJSETRNG Returns the original stream and sets the new RNG                                                                 
%                                                                        
% Usage: [] = vjSetRNG(rngStream)                                          
%                                                                        
% Expects:                                                               
% ----------------------------------------                               
% 'rng'       : Set default stream to this value                                                     
% 'seed'      : Set RNG using this seed 
% 
% Returns:                                                               
% ----------------------------------------                               
% orig        : The original RNG stream
%                                                                        
% Copyright (c) , Vinay Jethava (vjethava@gmail.com)                    
%                                                                        
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
[seed, varargin] = process_options(varargin, 'seed', 0); 
[rng, varargin] = process_options(varargin, 'rng', []); 
prevStream = RandStream.getDefaultStream(); 
if ~isempty(rng)
  RandStream.setDefaultStream(rng);
  fprintf(2, 'setRNG(): Reset stream to \n');
  rng.display();
elseif seed > 0
  s = RandStream.create('mt19937ar', 'seed', seed); 
  RandStream.setDefaultStream(s);
  fprintf(2, 'setRNG(): New stream with seed %d\n', seed); 
end
