function spm_standalone(varargin)
% A function to be compiled, which will run a standalone SPM.
%
% See http://www.mathworks.com/products/compiler/
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_standalone.m 4017 2010-07-26 17:14:31Z guillaume $ 

[v,r] = spm('Ver');
fprintf('%s (%s): %s\n',v,r,spm('Dir'));

if nargin && strcmpi(varargin{1},'run')
    spm('asciiwelcome');
    %spm('defaults','fmri');
    spm_jobman('initcfg');
    if nargin == 1
        h = spm_jobman;
        waitfor(h,'Visible','off');
    else
        for i=2:nargin
            try
                spm_jobman('run',varargin{i});
            catch
                fprintf('Execution of batch file ''%s'' failed.\n',varargin{i});
            end
        end
    end
    spm('Quit');
else
    spm(varargin{:});
end
