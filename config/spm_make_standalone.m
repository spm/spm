function spm_make_standalone(outdir, gateway, contentsver)
% Compile SPM as a standalone executable using the MATLAB Compiler
%   https://www.mathworks.com/products/compiler.html
%
% This will generate a standalone application, which can be run outside
% MATLAB, and therefore does not require a MATLAB licence.
%
% On Windows:
%   spm12.exe <modality>
%   spm12.exe batch <batch.m(at)>
%
% On Linux/Mac:
%   ./run_spm12.sh <MCRroot> <modality>
%   ./run_spm12.sh <MCRroot> batch <batch.m(at)>
%
% The first command starts SPM in interactive mode with GUI. The second
% executes a batch file or starts the Batch Editor if none is provided.
%
% Full list of options is accessible from:
%   ./run_spm12.sh <MCRroot> --help
%
% When deployed, compiled applications will require the MATLAB Runtime:
%   https://www.mathworks.com/products/compiler/matlab-runtime.html
% 
% See spm_standalone.m and https://en.wikibooks.org/wiki/SPM/Standalone
%__________________________________________________________________________
% Copyright (C) 2010-2019 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_make_standalone.m 7534 2019-02-20 17:09:45Z guillaume $


%-Check startup.m
%--------------------------------------------------------------------------
if exist('startup','file')
    warning('A startup.m has been detected in %s.\n',...
        fileparts(which('startup')));
end

%-Input arguments
%--------------------------------------------------------------------------
if ~nargin
    outdir = fullfile(spm('Dir'),'..','standalone'); 
    if ~exist(outdir,'dir'), mkdir(outdir); end
end
if nargin < 2 || isempty(gateway), gateway = 'spm_standalone.m'; end
if nargin < 3, contentsver = ''; end

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
        [dummy,I] = sort({f(9:end),version('-release')});
        nrmv = I(1) == 2;
    end
    if ~nrmv
        [sts, msg] = rmdir(d{i},'s');
    end
end

%==========================================================================
%-Compilation
%==========================================================================
Nopts = {'-p',fullfile(matlabroot,'toolbox','signal')};
if ~exist(Nopts{2},'dir'), Nopts = {}; end
Ropts = {'-R','-singleCompThread'} ;
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
