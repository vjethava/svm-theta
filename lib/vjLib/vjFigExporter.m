function [f] = vjFigExporter(fig_handle, varargin)
%% VJFIGEXPORTER This function exports figure in two formats (eps and fig)  
%   into the default directory path.  
%  
%   Options:
%   --------
%   'expName'     experiment name (default: calling function)
%   'project'     project name (used in identifying output directory)
%   'figDir'      output figure directory name
%   'resDir'      absolute directory path for figure directory (default:
%                   chosen based on project name)
%   'figNum'      figure number (in experiment) (default: 1) 
%   'log_handle'     log file handle (default: 2 (stderr) )
%   fig_handle    figure handle to save
%
%   @Note PDF export from matlab is not tight (exports to a4 page). Better
%   to export to epsc and then convert.
%
%   (c) Vinay Jethava, 2011
%   Chalmer University of Technology
%   
[project, varargin]  = process_options(varargin, 'project', 'InferGN');
[figDir, varargin] = process_options(varargin, 'figDir', 'figures'); 
currDir = cd ; %current Directory
baseDirPos = strfind(currDir, project); % find project name in full path
if size(baseDirPos) == 0 % string not found
  outDir = [currDir filesep];
else
  outDir0 = currDir(1:(baseDirPos(1)-1)); 
  outDir = [outDir0 filesep project filesep]; 
end
[resDir, varargin] = process_options(varargin, 'resDir', outDir); 
myStack  = dbstack; 
if (max(size(myStack)) > 1) & (isfield(myStack(2), 'name'))
  % use experiment name from calling function
  [expName] = process_options(varargin, 'expName', myStack(2).name); 
else
  % use randomly chosen file name
  [pathstr, tmp_name, ext, versn] = fileparts(tempname); 
  [expName] = process_options(varargin, 'expName', tmp_name); 
end
[figNum, varargin] = process_options(varargin, 'figNum', 1); 
[log_handle, varargin] = process_options(varargin, 'log_handle', 2); 
%% start the output
if exist([resDir filesep figDir],'dir') ~= 7
  mkdir(resDir, figDir); 
end
fprintf(log_handle, 'vjFigExporter() expName: %s figNum: %d\n', expName, figNum); 
% [outDir, varargin] = process_options(varargin, 'outDir', [resDir filesep figDir]);
fileName = [outDir filesep expName 'Fig' num2str(figNum)]; 
saveas(fig_handle, fileName, 'epsc'); 
saveas(fig_handle, fileName, 'fig');  
