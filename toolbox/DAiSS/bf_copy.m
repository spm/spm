function out = bf_copy
% Sets up a new analysis by copying an existing one
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_copy.m 7703 2019-11-22 12:06:29Z guillaume $

% dir Directory
% ---------------------------------------------------------------------

BF = cfg_files;
BF.tag = 'BF';
BF.name = 'BF.mat file';
BF.filter = '^BF.mat$';
BF.num = [1 1];
BF.help = {'Select BF.mat file.'};

dir         = cfg_files;
dir.tag     = 'dir';
dir.name    = 'Output directory';
dir.help    = {'Select a directory where the copied BF.mat file will be written.'};
dir.filter = 'dir';
dir.ufilter = '.*';
dir.num     = [1 1];

steps = cfg_menu;
steps.tag = 'steps';
steps.name = 'Steps to copy';
steps.help = {'Select ''all'' or the last step to copy'};
steps.labels = ['all'; bf_std_fields];
steps.values = ['all'; bf_std_fields];
steps.val = {'all'};

out = cfg_exbranch;
out.tag = 'copy';
out.name = 'Copy analysis';
out.val = {BF, dir, steps};
out.help = {'Make a copy of existing anallysis'};
out.prog = @bf_copy_run;
out.vout = @bf_copy_vout;
out.modality = {'EEG'};
end

function  out = bf_copy_run(job)

BF     = job.BF{1};
outdir = job.dir{1};

cd(outdir);

%-Ask about overwriting files from previous analyses
%--------------------------------------------------------------------------
if exist(fullfile(pwd,'BF.mat'),'file')
    str = {'Current directory contains existing BF file:',...
        'Continuing will overwrite existing file!'};
    if spm_input(str,1,'bd','stop|continue',[1,0],1,mfilename);
        fprintf('%-40s: %30s\n\n',...
            'Abort...   (existing BF file)',spm('time'));
        out = []; return
    end
end

if isequal(job.steps, 'all')
    [r, msg] = copyfile(BF, ...
        fullfile(outdir, 'BF.mat'), 'f');
    if ~r
        error(msg);
    end
else
    fields = bf_std_fields;
    ind    = strmatch(job.steps, fields);
    BF     = bf_load(BF, fields(1:ind));
    bf_save(BF, 1);
end

out.BF{1} = fullfile(outdir, 'BF.mat');

end

function dep = bf_copy_vout(job)
% Output is always in field "BF", no matter how job is structured
dep = cfg_dep;
dep.sname = 'BF.mat file';
% reference field "B" from output
dep.src_output = substruct('.','BF');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end
