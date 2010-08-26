function spm_print(varargin)
% Print the Graphics window
% FORMAT spm_print
% Print the default window to the default output file.
% FORMAT spm_print(filename)
% Print the default window to the specified output file.
% FORMAT spm_print('',figurename)
% Print the window with the specified name to the default output file.
% FORMAT spm_print(filename,figurename)
% Print the window with the specified name to the specified output file.
% FORMAT spm_print('',figurehandle)
% Print the figure with the specified handle to the default output file.
% FORMAT spm_print(filename,figurehandle)
% Print the figure with the specified handle to the specified output file.
% FORMAT spm_print(job)
% Run a batch print job, print the specified window to the specified
% output file.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_print.m 4050 2010-08-26 18:09:18Z guillaume $

%-Run spm_print through the Batch System to get configured print options
%==========================================================================
if ~nargin
    spm_jobman('serial','','spm.util.print','',{2},NaN);
    return;
elseif ischar(varargin{1})
    if nargin == 1
        spm_jobman('serial','','spm.util.print',varargin{1},{2},NaN);
    else
        spm_jobman('serial','','spm.util.print',varargin{1},{2}, ...
                   varargin{2});
    end
    return;
end

%-Print the Graphics window
%==========================================================================
job = varargin{1};
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
    if isfield(job.fig,'fighandle')
        if ishandle(job.fig.fighandle) && ...
                strcmpi(get(job.fig.fighandle,'type'),'figure')
            fg = job.fig.fighandle;
        elseif ~isfinite(job.fig.fighandle)
            fg = get(0,'CurrentFigure');
            if ~strcmp(get(fg,'Tag'),'Help')
                fg = spm_figure('FindWin','Graphics');
            end
        else
            fg = [];
        end
    else
        fg = spm_figure('FindWin',job.fig.figname);
    end
    if isempty(fg)
        fprintf('\nFigure not found: nothing has been printed.\n');
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
        fprintf('\nPrinting to\n%s%s%s\n',pwd,filesep,nam1);
    else
        fprintf('\nPrinting to\n%s\n',nam1);
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
