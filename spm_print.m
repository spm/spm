function spm_print(job)
% Print the graphics window
%____________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_print.m 2283 2008-10-01 14:25:09Z john $

% Run spm_print always as job to get configured print options
if nargin == 0
    spm_jobman('serial','','spm.util.print','');
    return;
elseif ischar(job)
    spm_jobman('serial','','spm.util.print',job);
    return;
end;


try
    mon = {'Jan','Feb','Mar','Apr','May','Jun',...
            'Jul','Aug','Sep','Oct','Nov','Dec'};
    t   = clock;
    nam = ['spm_' num2str(t(1)) mon{t(2)} sprintf('%.2d',t(3))];

    if isempty(job.fname)
        if job.opts.append,
            nam1 = fullfile(pwd,[nam job.opts.ext]);
        else
            nam1 = sprintf('%s_%3d',nam,1);
            for i=1:100000,
                nam1 = fullfile(pwd,sprintf('%s_%.3d%s',nam,i,job.opts.ext));
                if ~exist(nam1,'file'), break; end;
            end;
        end;
    else
        nam1 = job.fname;
    end;
    opts = {nam1,'-noui','-painters',job.opts.opt{:}};
    if strcmp(get(gcf,'Tag'),'Help'),
        fg = gcf;
    else
        fg = spm_figure('FindWin','Graphics');
    end;
    print(fg,opts{:});
    if isempty(strfind(nam1,filesep))
    fprintf('\nPrinting Graphics Windows to\n%s%s%s\n',pwd,filesep,nam1);
    else
    fprintf('\nPrinting Graphics Windows to\n%s\n',nam1);
    end
catch
    errstr = lasterr;
    tmp = [find(abs(errstr)==10),length(errstr)+1];
    str = {errstr(1:tmp(1)-1)};
    for i = 1:length(tmp)-1
        if tmp(i)+1 < tmp(i+1)
            str = [str, {errstr(tmp(i)+1:tmp(i+1)-1)}];
        end
    end
    str = {str{:},  '','- Print options are:', opts{:},...
                    '','- Current directory is:',['    ',pwd],...
                    '','            * nothing has been printed *'};
    spm('alert!',str,'printing problem...',sqrt(-1));
end;

