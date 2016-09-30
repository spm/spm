function varargout = spm_print(varargin)
% Print figure 
% FORMAT spm_print(fname,F,opts)
% fname  - output filename [Default: 'spm_<date>']
% F      - figure handle or tag [Default: 'Graphics']
% opts   - structure containing printing options
%          [Default: defaults.ui.print from spm_defaults.m]
%
% FORMAT spm_print(job)
% Run a batch print job (see spm_cfg_print)
%__________________________________________________________________________
% Copyright (C) 1994-2016 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_print.m 6894 2016-09-30 16:48:46Z spm $


%-Shortcut for list of graphics file formats available
%--------------------------------------------------------------------------
if nargin > 0 && isequal(lower(varargin{1}),'format')
    [varargout{1:nargout}] = print_format(varargin{2:end});
    return;
end

%-Conversion layer for spm.util.print batch
%--------------------------------------------------------------------------
if nargin == 1 && isstruct(varargin{1})
    job = varargin{1};
    spm_print(job.fname, job.fig.(char(fieldnames(job.fig))), job.opts);
    return;
end

%-Input parameters
%--------------------------------------------------------------------------
if nargin<1, fname = ''; else    fname = varargin{1}; end

if nargin<2, F = NaN;    else    F     = varargin{2}; end
if isnumeric(F) && ~isfinite(F), F     = [];          end
if isempty(F),                   F     = 'Graphics';  end

opts = spm_get_defaults('ui.print');
if ischar(opts)
    opts = spm_print('format',opts);
    if isempty(opts), disp('Print warning: unknown printing format.'); end
end
if nargin == 3
    o = varargin{3};
    if ischar(o), o = spm_print('format',o); end
    if isempty(o), disp('Print warning: unknown printing format.'); end
    try, opts.ext    = o.ext;    end
    try, opts.append = o.append; end
    try, opts.opt    = o.opt;    end
end
opts.opt = [opts.opt {'-noui','-painters'}];
if opts.append, opts.opt = [opts.opt {'-append'}]; end

%-Get output filename
%--------------------------------------------------------------------------
if isempty(fname)
    bname = ['spm_' datestr(now,'yyyymmmdd')];
    pname = pwd;
else
    bname = spm_file(fname,'basename');
    pname = spm_file(fname,'fpath');
end
fname = fullfile(pname,[bname opts.ext]);
if ~opts.append
    fname = spm_file(fname,'unique');
end

%-Get figure handle
%--------------------------------------------------------------------------
F = spm_figure('FindWin',F);
%if isempty(F), F = get(0,'CurrentFigure'); end
if isempty(F)
    disp('Print error: Figure not found.');
    return;
end

%-Print
%==========================================================================
try
    if ismember('-dfig',opts.opt)
        saveas(F, fname, 'fig');
    else
        if isdeployed
            deployprint(F, fname, opts.opt{:});
        else
            print(F, fname, opts.opt{:});
        end
    end
catch
    disp('Print error: nothing has been printed.');
    l = lasterror;
    disp(l.message);
    return;
end

%-Report
%--------------------------------------------------------------------------
str = '';
if ~isempty(get(F,'Tag')), str = sprintf('''%s'' ',get(F,'Tag')); end
fprintf('Printing %sfigure to:\n',str);
if ispc && ~ismember('-dfig',opts.opt), cmd = 'winopen(''%s'')';
else                                    cmd = 'open(''%s'')'; end
fprintf('  %s\n',spm_file(fname,'link',cmd));

if nargout, varargout = { fname }; end


%==========================================================================
% function pf = print_format(f)
%==========================================================================
function [pf, i] = print_format(f)

pf(1).name   = 'PostScript (PS)';
pf(1).opt    = {'-dpsc2'};
pf(1).append = true;
pf(1).ext    = '.ps';
pf(1).label  = {'ps'};

pf(2).name   = 'Encapsulated PostScript (EPS)';
pf(2).opt    = {'-depsc2'};
pf(2).append = false;
pf(2).ext    = '.eps';
pf(2).label  = {'eps'};

pf(3).name   = 'Portable Document Format (PDF)';
pf(3).opt    = {'-dpdf'};
pf(3).append = false;
pf(3).ext    = '.pdf';
pf(3).label  = {'pdf'};

pf(4).name   = 'JPEG image';
pf(4).opt    = {'-djpeg'};
pf(4).append = false;
pf(4).ext    = '.jpg';
pf(4).label  = {'jpg','jpeg'};

pf(5).name   = 'PNG image';
pf(5).opt    = {'-dpng'};
pf(5).append = false;
pf(5).ext    = '.png';
pf(5).label  = {'png'};

pf(6).name   = 'TIFF image';
pf(6).opt    = {'-dtiff'};
pf(6).append = false;
pf(6).ext    = '.tif';
pf(6).label  = {'tif','tiff'};

pf(7).name   = 'MATLAB figure';
pf(7).opt    = {'-dfig'}; % use saveas instead of print
pf(7).append = false;
pf(7).ext    = '.fig';
pf(7).label  = {'fig'};

if nargin
    for i=1:numel(pf)
        if ismember(lower(f),lower(pf(i).label))
            pf = pf(i); return;
        end
    end
    pf = []; i = [];
end
