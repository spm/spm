function spm_make_standalone(outdir, gateway, contentsver, tbxs)
% Compile SPM as a standalone executable using the MATLAB Compiler
%   https://www.mathworks.com/products/compiler.html
%
% This will generate a standalone application, which can be run outside
% MATLAB, and therefore does not require a MATLAB licence.
%
% On Windows:
%   spm.exe <modality>
%   spm.exe batch <batch.m(at)>
%   spm.exe script script.m
%
% On Linux/Mac:
%   ./run_spm.sh <MCRroot> <modality>
%   ./run_spm.sh <MCRroot> batch <batch.m(at)>
%   ./run_spm.sh <MCRroot> script script.m
%
% The first command starts SPM in interactive mode with GUI. The second
% executes a batch file or starts the Batch Editor if none is provided,
% while the third command evaluates the content of script.m. Extra
% command line arguments are available in a cell array variable named
% "inputs".
%
% Full list of options is accessible from:
%   ./run_spm.sh <MCRroot> --help
%
% When deployed, compiled applications will require the MATLAB Runtime:
%   https://www.mathworks.com/products/compiler/matlab-runtime.html
% 
% See https://www.fil.ion.ucl.ac.uk/spm/docs/installation/standalone/ and
% spm_standalone.m
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2010-2023 Wellcome Centre for Human Neuroimaging


%-Check startup.m
%--------------------------------------------------------------------------
if exist('startup','file')
    warning('A startup.m has been detected in %s.\n',...
        fileparts(which('startup')));
end

%-Input arguments
%--------------------------------------------------------------------------
if ~nargin || isempty(outdir)
    outdir = fullfile(spm('Dir'),'..','standalone'); 
end
if ~exist(outdir,'dir'), mkdir(outdir); end
if nargin < 2 || isempty(gateway), gateway = 'spm_standalone.m'; end
if nargin < 3, contentsver = ''; end
if nargin < 4, tbxs = {'signal'}; end

%==========================================================================
%-Static listing of SPM toolboxes
%==========================================================================
fid = fopen(fullfile(spm('Dir'),'config','spm_cfg_static_tools.m'),'wt');
fprintf(fid,'function values = spm_cfg_static_tools\n');
fprintf(fid,...
    '%% Static listing of all batch configuration files in the SPM toolbox folder\n');
%-Get the list of toolbox directories
tbxdir = fullfile(spm('Dir'),'toolbox');
d = [tbxdir; cellstr(spm_select('FPList',tbxdir,'dir'))];
ft = {};
%-Look for '*_cfg_*.m' files in these directories
for i=1:numel(d)
    fi = spm_select('List',d{i},'.*_cfg_.*\.m$');
    if ~isempty(fi)
        ft = [ft(:); cellstr(fi)];
    end
end
%-Create code to insert toolbox config
if isempty(ft)
    ftstr = '';
else
    ft = spm_file(ft,'basename');
    ftstr = sprintf('%s ', ft{:});
end
fprintf(fid,'values = {%s};\n', ftstr);
fclose(fid);

%==========================================================================
%-Static listing of batch application initialisation files
%==========================================================================
cfg_util('dumpcfg');

%==========================================================================
%-Duplicate Contents.m in Contents.txt for use in spm('Ver')
%==========================================================================
sts = copyfile(fullfile(spm('Dir'),'Contents.m'),...
               fullfile(spm('Dir'),'Contents.txt'));
if ~sts, warning('Copy of Contents.m failed.'); end
if ~isempty(contentsver)
    % Format: 'xxxx (SPMx) dd-mmm-yyyy'
    f = fileread(fullfile(spm('Dir'),'Contents.txt'));
    f = regexprep(f,'% Version \S+ \S+ \S+',['% Version ' contentsver]);
    fid = fopen(fullfile(spm('Dir'),'Contents.txt'),'w');
    fprintf(fid,'%s',f);
    fclose(fid);
end

%==========================================================================
%-Trim FieldTrip
%==========================================================================
d = fullfile(spm('Dir'),'external','fieldtrip','compat');
d = cellstr(spm_select('FPList',d,'dir'));
for i=1:numel(d)
    f = spm_file(d{i},'basename');
    nrmv = strncmp(f,'matlablt',8);
    if nrmv
        [~,I] = sort({f(9:end),version('-release')});
        nrmv = I(1) == 2;
    end
    if ~nrmv
        [~] = rmdir(d{i},'s');
    end
end
for i=1:numel(tbxs)
    d = fullfile(spm('Dir'),'external','fieldtrip','external',tbxs{i});
    if exist(d,'dir')
        [~] = rmdir(d,'s');
    end
end

%==========================================================================
%-Compilation
%==========================================================================
Nopts = {};
for i=1:numel(tbxs)
    d = fullfile(matlabroot,'toolbox',tbxs{i});
    if exist(d,'dir')
        Nopts = [Nopts {'-p'} {d}];
    end
end
Ropts = {} ;
if ~ismac && spm_check_version('matlab','8.4') >= 0
    Ropts = [Ropts, {'-R','-softwareopengl'}];
end
mcc('-m', '-C', '-v',...
    '-o',lower(spm('Ver')),...
    '-d',outdir,...
    '-N',Nopts{:},...
    Ropts{:},...
    '-a',spm('Dir'),...
    gateway);
