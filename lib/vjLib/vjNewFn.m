function [] = vjNewFn(filename, inputs, outputs)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% VJNEWFN creates new function file with right structure                 %
%                                                                        %
% Usage: [outputs] = vjNewFn(filename, inputs, outputs)                  %
%                                                                        %
% Expects:                                                               %
% ----------------------------------------                               %
% filename          : Name of the function                               %
% inputs            : Cells containing list of inputs                    %
% outputs           : Cells containing list of outputs                   %
%                                                                        %
% Returns:                                                               %
% ----------------------------------------                               %
%                                                                        %
% Copyright (c) 2010, Vinay Jethava (vjethava@gmail.com)                 %
%                                                                        %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
if nargin < 2
  inputs = {'vargargin'};
  outputs = {};
end
if ~strcmpi(filename((end-1):end), '.m')
    filename = [filename '.m'];
end
if exist(filename, 'file')
    display('File Already Exists');
else
    % f = 2; %% write to stderr instead
    fnName = filename(1:(end-2)); 
    c = {}; % the comments
    % function title
    c{length(c) + 1} = ['function ' getHeader(fnName, inputs, outputs)]; 
    % Descriptive one-liner
    c{length(c) + 1} = getComment(getLine()); 
    c{length(c) + 1} = getComment(upper(fnName)); % Description line
    c{length(c) + 1} = getComment(); 
    % Usage example
    c{length(c) + 1} = getComment(['Usage: ' getHeader(fnName, inputs, outputs)]);
    c{length(c) + 1} = getComment(); 
    % Insert input block 
    d = getDescription('Expects:', inputs);
    for i=1:length(d)
      c{length(c) + 1} = getComment(d{i});
    end
    c{length(c) + 1} = getComment(); 
    % Insert output block
    d = getDescription('Returns:', outputs);
    for i=1:length(d)
      c{length(c) + 1} = getComment(d{i});
    end
    c{length(c) + 1} = getComment();
    % Insert Note block
    c{length(c) + 1} = getComment('Note:');
    c{length(c) + 1} = getComment(getLine('-', 40));
    c{length(c) + 1} = getComment();
    
    % Insert See Also block
    c{length(c) + 1} = getComment('See Also:');
    c{length(c) + 1} = getComment(getLine('-', 40));
    c{length(c) + 1} = getComment();
    
    % Insert Copyright block 
    c{length(c) + 1} = getComment(getCopyright());
    c{length(c) + 1} = getComment(); 
    c{length(c) + 1} = getComment(getLine());
    try 
        f = fopen(filename, 'w');
        for i=1:length(c)
          fprintf(f, '%s\n', c{i});
        end
    catch
        fclose(f);
    end
end
try
     open(filename);
catch
     display('Could not open created file!');
end
end %%% Main function ends here

function [line] = getLine(c, num)
  if nargin < 1
    c = '%';
    num = 70; 
  end
    line = sprintf('%s', char( c * ones(1, num)));
    % underLine = repmat(sprintf(uChar), 1, uLen); 
end

function [header]=getHeader(fnName, inputs, outputs)
%% GETHEADERNAME Returns the header line for the function
header = '[';
if (nargin >= 3) && (~isempty(outputs))
  so = length(outputs);
  for i = 1:(so-1)
    header = [ header outputs{i} ', '];
  end
  header = [header outputs{so}];
end
header = [header '] = ' fnName '('];
if (nargin >= 2) && (~isempty(inputs))
  si = length(inputs);
  for i = 1:(si-1)
    header = [ header inputs{i} ', '];
  end
  header = [header inputs{si}];
end
header = [header ')'];
end

function [comment] = getComment(line) 
  if nargin < 1
    line = ''; % empty line
  end
  comment = sprintf('%s%-70s%s', getPrefix(), line, getPostfix()); 
end 

function [copyright] = getCopyright(name, cyear, email)
  if nargin < 1
    name = 'Vinay Jethava';
    cyear = year(date);
    email= 'vjethava@gmail.com';
  end
  copyright = ['Copyright (c) ' sprintf('%d', cyear) ', ' name ' (' email ')'];
end

function [prefix] = getPrefix() 
  prefix = '% '; 
end
function [postfix] = getPostfix() 
  postfix = ' '; 
end


function [d] = getDescription(paramType, params, uChar)
%% GETDESCRIPTION Returns the description for the parameter
d = {};
if nargin < 3
  uChar = '-';
end
underLine = getLine(uChar, 40);
d{length(d) + 1} =  (paramType); 
d{length(d) + 1} =  underLine; 
if (nargin >= 2) && (~isempty(params))
  for i=1:length(params)
    d{length(d) + 1} =  sprintf('%-15s : ', params{i});  
  end
end
end