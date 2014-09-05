function spm_standalone(varargin)
% A function to be compiled, which will run a standalone SPM.
%
% See MATLAB Compiler: http://www.mathworks.com/products/compiler/
% See also config/spm_make_standalone.m
%__________________________________________________________________________
% Copyright (C) 2010-2014 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_standalone.m 6156 2014-09-05 17:34:53Z guillaume $ 

[v,r] = spm('Ver');
fprintf('%s (%s): %s\n',v,r,spm('Dir'));

if ~nargin, action = ''; else action = varargin{1}; end

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
        if nargin > 1
            spm('Run',varargin(2));
        else
            spm('Run');
        end
        
    otherwise
    %----------------------------------------------------------------------
        spm(varargin{:});
        
end
