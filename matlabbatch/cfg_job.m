classdef cfg_job
% A command line frontend to the batch system
% obj = cfg_job
% The cfg_job object allows to create batch jobs interactively on the
% MATLAB command line. It will allow any cell/struct assignments and
% references that are allowed by the current batch configuration.
% As an option, for each subscript reference hints about eligible values
% can be displayed. Note that hints display disables MATLABs
% autocompletion feature.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_job.m 6460 2015-05-28 08:30:28Z volkmar $

    properties (Access=private)
        c0;
        showhints = false;
    end
    methods
        %% Constructor, Access to internal properties
        function obj = cfg_job(c0)
        % Initialise object
        % obj = cfg_job
        % Initialise object using information from cfg_util.
        % obj = cfg_job(c0)
        % Initialise object using an explicitly specified cfg_item
        % configuration tree.
            if nargin == 0
                obj.c0 = cfg_util('getcfg');
            else
                obj.c0 = c0;
            end
        end
        function obj = hints(obj, val)
        % Change hint display settings
        % hints(obj, flag)
        % If flag is true, hints are displayed each time a subscript
        % reference is run. These hints include information about
        % eligible values, current values etc for the referenced
        % cfg_item.
        % Note that hints display interferes with MATLAB autocompletion
        % mechanism.
            if nargin == 1
                obj.showhints = ~obj.showhints;
            else
                if islogical(val)
                    obj.showhints = val;
                elseif ischar(val)
                    switch lower(val)
                        case {'on','yes','true'}
                            obj.showhints = true;
                        case {'off','no','false'}
                            obj.showhints = false;
                    end
                end
            end
        end
        %% Overridden subscript assignment and reference
        % The core idea is to intercept MATLABs subsasgn/subsref
        % calls. These are redirected to cfg_item methods that traverse a
        % configuration tree based on job subscripts.
        function obj = subsasgn(obj, subs, val)
            obj.c0 = subsasgn_job(obj.c0, subs, val);
        end
        function varargout = subsref(obj, subs)
            [ritem, varargout{1:max(nargout,1)}] = subsref_job(obj.c0, subs, obj.c0);
            if obj.showhints || nargout == 0
                str = showdetail(ritem);
                fprintf('%s\n', str{:});
            end
        end
        %% Disp method
        function disp(obj)
            str = showdetail(obj.c0);
            fprintf('%s\n', str{:});
            [~, val] = harvest(obj.c0, obj.c0, false, false);
            fprintf('\nCurrent contents:\n')
            disp(val)
        end
        %% Job specific methods
        function [tag, val] = harvest(obj)
        % Harvest a job.
        % [tag, val] = harvest(job)
        % This returns the tag of the top level cfg_item object and a
        % MATLAB struct/cell variable containing the job inputs. 
            [tag, val] = harvest(obj.c0, obj.c0, false, false);
        end
        function obj = initialise(obj, val)
        % Initialise obj with a job.
        % initialise(obj, job)
        % Load a MATLAB struct/cell job into a cfg_obj.
            obj.c0 = initialise(obj.c0, '<DEFAULTS>', false);
            obj.c0 = initialise(obj.c0, val, false);
        end
        function savejob(obj, filename)
        % Save job variable to .m file.
        % savejob(obj, filename)
        % Save the harvested job to a MATLAB script.
            [~, matlabbatch] = harvest(obj);
            str = gencode(matlabbatch);
            [fid, msg] = fopen(filename, 'w');
            if fid == -1
                cfg_message('matlabbatch:fopen', 'Failed to open ''%s'' for writing:\n%s', filename, msg);
            end
            fprintf(fid, 'matlabbatch = %s;\n', class(obj));
            fprintf(fid, '%s\n', str{:});
            fclose(fid);
        end
    end
end
