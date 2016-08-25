function spm_standalone(varargin)
% A function to be compiled, which will run a standalone SPM.
%
% See MATLAB Compiler: http://www.mathworks.com/products/compiler/
% See also config/spm_make_standalone.m
%__________________________________________________________________________
% Copyright (C) 2010-2016 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_standalone.m 6862 2016-08-25 14:42:19Z guillaume $ 

[v,r] = spm('Ver');
fprintf('%s (%s): %s\n',v,r,spm('Dir'));

if ~nargin, action = ''; else action = varargin{1}; end

exit_code = 0;

switch lower(action)
    
    case {'batch', 'run'}
    %----------------------------------------------------------------------
        spm('asciiwelcome');
        %spm('defaults','fmri');
        spm_jobman('initcfg');
        if nargin == 1
            h = spm_jobman;
            waitfor(h,'Visible','off');
        else
            %spm_get_defaults('cmdline',true);
            for i=2:nargin
                try
                    spm_jobman('run',varargin{i});
                catch
                    fprintf('Execution failed: %s', varargin{i});
                end
            end
        end
        spm('Quit');
        
    case 'script'
    %----------------------------------------------------------------------
        spm('asciiwelcome');
        assignin('base','inputs',varargin(3:end));
        try
            if nargin > 1
                spm('Run',varargin(2));
            else
                spm('Run');
            end
        catch
            exit_code = 1;
        end
        
    case 'function'
    %----------------------------------------------------------------------
        spm('asciiwelcome');
        if nargin == 1
            fcn = spm_input('function name','!+1','s','');
        else
            fcn = varargin{2};
        end
        try
            feval(fcn,varargin{3:end});
        catch
            exit_code = 1;
        end
    
    otherwise
    %----------------------------------------------------------------------
        spm(varargin{:});
        
end

%-Display error message and return exit code (or use finish.m script)
%--------------------------------------------------------------------------
if exit_code ~= 0
    err = lasterror;
    msg{1} = err.message;
    if isfield(err,'stack')
        for i=1:numel(err.stack)
            if err.stack(i).line
                l = sprintf(' (line %d)',err.stack(i).line);
            else
                l = '';
            end
            msg{end+1} = sprintf('Error in %s%s',err.stack(i).name,l);
        end
    end
    fprintf('%s\n',msg{:});
    
    exit(exit_code);
end
