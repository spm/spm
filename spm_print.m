function spm_print(job)
% Print the Graphics window
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_print.m 3933 2010-06-17 14:19:08Z guillaume $

%-Run spm_print through the Batch System to get configured print options
%==========================================================================
if ~nargin
    spm_jobman('serial','','spm.util.print','');
    return;
elseif ischar(job)
    spm_jobman('serial','','spm.util.print',job);
    return;
end

%-Print the Graphics window
%==========================================================================
try
    %-Get output filename
    %----------------------------------------------------------------------
    if isempty(job.fname)
        nam = ['spm_' datestr(now,'yyyymmmdd')];
        if job.opts.append
            nam1 = fullfile(pwd,[nam job.opts.ext]);
        else
            nam1=''; i=1;
            while isempty(nam1)
                nam1 = fullfile(pwd,sprintf('%s_%.3d%s',nam,i,job.opts.ext));
                if ~exist(nam1,'file'), break; else nam1='';i=i+1; end;
            end
        end
    else
        nam1 = job.fname;
    end
    
    %-Get print options
    %----------------------------------------------------------------------
    opts = {nam1,'-noui','-painters',job.opts.opt{:}};
    
    %-Get figure handle
    %----------------------------------------------------------------------
    if strcmp(get(gcf,'Tag'),'Help')
        fg = gcf;
    else
        fg = spm_figure('FindWin','Graphics');
    end
    if isempty(fg)
        fprintf('\nGraphics window not found: nothing has been printed.\n');
        return;
    end
    
    %-Print
    %----------------------------------------------------------------------
    if isdeployed
        deployprint(fg,opts{:});
    else
        print(fg,opts{:});
    end
    
    %-Report
    %----------------------------------------------------------------------
    if isempty(strfind(nam1,filesep))
        fprintf('\nPrinting Graphics Windows to\n%s%s%s\n',pwd,filesep,nam1);
    else
        fprintf('\nPrinting Graphics Windows to\n%s\n',nam1);
    end
    
%-Report errors
%==========================================================================
catch
    errstr = lasterror;
    errstr = errstr.message;
    str    = textscan(errstr,'%s','delimiter',sprintf('\n'));
    str    = str{1};
    str    = {str{:},'','- Print options are:', opts{:},...
                     '','- Current directory is:',['    ',pwd],...
                     '','            * nothing has been printed *'};
    spm('alert!',str,'printing problem...',sqrt(-1));
end
