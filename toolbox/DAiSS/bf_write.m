function out = bf_write
% Writes out the results of beamforming analysis
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_write.m 7703 2019-11-22 12:06:29Z guillaume $

% dir Directory
% ---------------------------------------------------------------------
BF = cfg_files;
BF.tag = 'BF';
BF.name = 'BF.mat file';
BF.filter = '^BF.mat$';
BF.num = [1 1];
BF.help = {'Select BF.mat file.'};

%--------------------------------------------------------------------------
% method
%--------------------------------------------------------------------------
plugin      = cfg_choice;
plugin.tag  = 'plugin';
plugin.name = 'What to write out';

write_funs = spm_select('List', fileparts(mfilename('fullpath')), '^bf_write_.*\.m$');
write_funs = cellstr(write_funs );
for i = 1:numel(write_funs)
    plugin.values{i} = feval(spm_file(write_funs{i},'basename'));
end

out = cfg_exbranch;
out.tag = 'write';
out.name = 'Write';
out.val = {BF, plugin};
out.help = {'Write out results'};
out.prog = @bf_write_run;
out.vout = @bf_write_vout;
out.modality = {'EEG'};
end

function  out = bf_write_run(job)

outdir = spm_file(job.BF{1}, 'fpath');

cd(outdir);

BF = bf_load('BF.mat', {'data', 'sources', 'features', 'output'});

plugin_name   = cell2mat(fieldnames(job.plugin));

outfield_name =  strtok(plugin_name, '_');

BF.write.(outfield_name) = feval(['bf_write_' plugin_name], BF, job.plugin.(plugin_name));

bf_save(BF);

out.BF{1} = fullfile(outdir, 'BF.mat');
out.files = BF.write.(outfield_name).files;
end

function dep = bf_write_vout(job)
% Output is always in field "D", no matter how job is structured
dep(1) = cfg_dep;
dep(1).sname = 'BF.mat file';
% reference field "B" from output
dep(1).src_output = substruct('.','BF');
% this can be entered into any evaluated input
dep(1).tgt_spec   = cfg_findspec({{'filter','mat'}});
dep(2).sname      = 'Output files';
dep(2).src_output = substruct('.','files');

plugin_name   = cell2mat(fieldnames(job.plugin));

switch plugin_name
    case 'spmeeg'
        dep(2).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
    case 'nifti'
        dep(2).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    case 'gifti'
        dep(2).tgt_spec   = cfg_findspec({{'filter','mesh','strtype','e'}});
end
end