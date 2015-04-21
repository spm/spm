function spm_make_standalone(outdir)
% Compile SPM as a standalone executable using the MATLAB compiler
%   http://www.mathworks.com/products/compiler/
%
% This will generate a standalone program, which can be run
% outside MATLAB, and therefore does not use up a MATLAB licence.
%
% On Windows:
%   spm12_wxx.exe <modality>
%   spm12_wxx.exe run <batch.m(at)>
%
% On Linux/Mac:
%   ./run_spm12.sh <MCRroot> <modality>
%   ./run_spm12.sh <MCRroot> run <batch.m(at)>
%
% The first command starts SPM in interactive mode with GUI. The second
% executes a batch file or starts the Batch Editor if empty.
%
% See spm_standalone.m
%__________________________________________________________________________
% Copyright (C) 2010-2015 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_make_standalone.m 6416 2015-04-21 15:34:10Z guillaume $

%-Care of startup.m
%--------------------------------------------------------------------------
% see http://www.mathworks.com/support/solutions/data/1-QXFMQ.html?1-QXFMQ
if exist('startup','file')
    warning('A startup.m has been detected in %s.\n',...
        fileparts(which('startup')));
end

if ~nargin, outdir = fullfile(spm('dir'),'..','spm_exec'); mkdir(outdir); end

%==========================================================================
%-Static listing of SPM toolboxes
%==========================================================================
fid = fopen(fullfile(spm('dir'),'config','spm_cfg_static_tools.m'),'wt');
fprintf(fid,'function values = spm_cfg_static_tools\n');
fprintf(fid,...
    '%% Static listing of all batch configuration files in the SPM toolbox folder\n');
% create code to insert toolbox config
%-Toolbox autodetection
%-Get the list of toolbox directories
tbxdir = fullfile(spm('Dir'),'toolbox');
d  = dir(tbxdir); d = {d([d.isdir]).name};
dd = regexp(d,'^\.');
%(Beware, regexp returns an array if input cell array is of dim 0 or 1)
if ~iscell(dd), dd = {dd}; end
d  = {'' d{cellfun('isempty',dd)}};
ft = {};
%-Look for '*_cfg_*.m' files in these directories
for i=1:length(d)
    d2 = fullfile(tbxdir,d{i});
    di = dir(d2); di = {di(~[di.isdir]).name};
    f2 = regexp(di,'.*_cfg_.*\.m$');
    if ~iscell(f2), f2 = {f2}; end
    fi = di(~cellfun('isempty',f2));
    if ~isempty(fi)
        ft = [ft(:); fi(:)];
    end
end
if ~isempty(ft)
    if isempty(ft)
        ftstr = '';
    else
        ft = cellfun(@(cft)strtok(cft,'.'),ft,'UniformOutput',false);
        ftstr  = sprintf('%s ', ft{:});
    end
    fprintf(fid,'values = {%s};\n', ftstr);
end
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

%==========================================================================
%-Compilation
%==========================================================================
opts = {'-p',fullfile(matlabroot,'toolbox','signal')};
mcc('-m', '-C', '-v',...
    '-o',lower(spm('Ver')),...
    '-d',outdir,...
    '-N',opts{:},...
    '-R','-singleCompThread',...
    '-a',spm('Dir'),...
    'spm_standalone.m');
