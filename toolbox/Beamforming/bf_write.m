function out = bf_write
% Writes out the results of beamforming analysis
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_write.m 4897 2012-09-04 16:32:18Z vladimir $

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

write_funs = spm_select('List', fullfile(spm('dir'), 'toolbox', 'Beamforming'), '^bf_write_.*\.m$');
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

BF = bf_load('BF.mat', {'data', 'sources', 'output'});

plugin_name   = cell2mat(fieldnames(job.plugin));

outfield_name =  strtok('plugin_name', '_');

BF.write.(outfield_name) = feval(['bf_write_' plugin_name], BF, job.plugin.(plugin_name));

bf_save(BF);

out.BF{1} = fullfile(outdir, 'BF.mat');
end

function dep = bf_write_vout(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'BF.mat file';
% reference field "B" from output
dep.src_output = substruct('.','BF');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});
end
