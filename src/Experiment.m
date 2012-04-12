classdef Experiment < handle
  properties 
    name
    stats
    caption
    precision=4; 
  end
  methods (Static)
    function [table]=get_table(stats)
      fields = fieldnames(stats);
      num_rows = length(stats); 
      num_fields = length(fieldnames(stats)); 
      table = [];
      for j=1:(num_rows+1)
        for i=1:num_fields
          if j==1
            table = [table sprintf('%s', fields{i})];
          else
            val = getfield(stats,{j-1}, fields{i});
            if ischar(val)
              table = [table sprintf('%s', val)]; 
            else
              table = [table sprintf('$%-.4g$', val)];
            end
          end
          if  (i < num_fields)
            table = [table sprintf('\t& ')];
          elseif i == num_fields
            table = [table sprintf(' \\\\\\hline\n')];
          end
        end        
      end
    end
  end
  properties (Dependent, SetAccess='private')
    num_rows
    num_fields
  end
  properties (Dependent)
    table
    figure
  end
  methods
    function [nr]=get.num_rows(obj)
      nr = length(obj.stats); 
    end
    function [nf]=get.num_fields(obj)
      nf = length(fieldnames(obj.stats)); 
    end
    function [table]=get.table(obj)
      table = Experiment.get_table(obj.stats); 
    end
    function disp(obj)
      if ~isempty(obj.name)
        fprintf(1, '\nExperiment: %s\n', obj.name);        
      end
      if ~isempty(obj.caption)
        fprintf(1, '%s\n', obj.caption);        
      end
      fprintf(1, '\n%s', obj.table); 
    end
  end
end
