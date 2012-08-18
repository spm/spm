function out = bf_postprocessing
% Performs postprocessing based on beamforming projectors
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_postprocessing.m 4849 2012-08-18 12:51:28Z vladimir $

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
plugin.name = 'Postprocessing method';

postprocessing_funs = spm_select('List', fullfile(spm('dir'), 'toolbox', 'Beamforming'), '^bf_postprocessing_.*\.m$');
postprocessing_funs = cellstr(postprocessing_funs );
for i = 1:numel(postprocessing_funs)
    plugin.values{i} = feval(spm_file(postprocessing_funs{i},'basename'));
end

out = cfg_exbranch;
out.tag = 'postprocessing';
out.name = 'Postprocessing';
out.val = {BF, plugin};
out.help = {'Perform postprocessing'};
out.prog = @bf_postprocessing_run;
out.vout = @bf_postprocessing_vout;
out.modality = {'EEG'};
end

function  out = bf_postprocessing_run(job)

outdir = spm_file(job.BF{1}, 'fpath');

cd(outdir);

BF = bf_load('BF.mat', {'data', 'inverse'});

plugin_name   = cell2mat(fieldnames(job.plugin));

outfield_name =  strtok(plugin_name, '_');

BF.postprocessing.(outfield_name) = feval(['bf_postprocessing_' plugin_name], BF, job.plugin.(plugin_name));

bf_save(BF);

out.BF{1} = fullfile(outdir, 'BF.mat');
end

function dep = bf_postprocessing_vout(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'BF.mat file';
% reference field "B" from output
dep.src_output = substruct('.','BF');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});
end
