function spm_standalone(varargin)
% A function to be compiled, which will run a standalone SPM.
%
% See MATLAB Compiler: http://www.mathworks.com/products/compiler/
% See also config/spm_make_standalone.m
%__________________________________________________________________________
% Copyright (C) 2010-2011 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_standalone.m 4237 2011-03-10 19:46:33Z guillaume $ 

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
        if nargin > 1
            inputs = varargin(3:end);
            fid = fopen(varargin{2});
            if fid == -1, error('Cannot open %s',varargin{2}); end
            S = fscanf(fid,'%c');
            fclose(fid);
            try
                eval(S);
            catch
                fprintf('Execution failed: %s\n',varargin{2});
                rethrow(lasterror);
            end
        end
        
    otherwise
    %----------------------------------------------------------------------
        spm(varargin{:});
        
end
