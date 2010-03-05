function exec_spm(varargin)
% A function to be compiled, which will run SPM.
%
% See http://www.mathworks.com/products/compiler/
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: exec_spm.m 3752 2010-03-05 12:54:50Z guillaume $ 

if nargin && strcmpi(varargin{1},'run')
    spm('asciiwelcome');
    spm_jobman('initcfg');
    if nargin == 1
        h = spm_jobman;
        waitfor(h,'Visible','off');
    else
        for i=2:nargin
            spm_jobman('run',varargin{i});
        end
    end
    spm('Quit');
else
    spm(varargin{:});
end
