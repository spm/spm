function out = bf_inverse
% Computes inverse projectors
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_inverse.m 7703 2019-11-22 12:06:29Z guillaume $

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
plugin.name = 'Inverse method';

inverse_funs = spm_select('List', fileparts(mfilename('fullpath')), '^bf_inverse_.*\.m$');
inverse_funs = cellstr(inverse_funs );
for i = 1:numel(inverse_funs)
    plugin.values{i} = feval(spm_file(inverse_funs{i},'basename'));
end


out = cfg_exbranch;
out.tag = 'inverse';
out.name = 'Inverse solution';
out.val = {BF, plugin};
out.help = {'Compute inverse projectors'};
out.prog = @bf_inverse_run;
out.vout = @bf_inverse_vout;
out.modality = {'EEG'};
end

function  out = bf_inverse_run(job)

outdir = spm_file(job.BF{1}, 'fpath');

cd(outdir);

BF = bf_load('BF.mat', {'data', 'sources', 'features'});
plugin_name = cell2mat(fieldnames(job.plugin));
S = job.plugin.(plugin_name);

if ~isa(S, 'struct')
    S = [];
end

modalities = intersect(fieldnames(BF.features), {'MEG', 'MEGPLANAR', 'MEGMAG', 'EEG'});

for m = 1:numel(modalities)
    S(1).modality = modalities{m};
    
    [S.L, channels] = bf_fuse_lf(BF, S.modality);
    
    BF.inverse.(modalities{m}) = feval(['bf_inverse_' plugin_name], BF, S);
    BF.inverse.(modalities{m}).channels = channels;
    if ~isfield(BF.inverse.(modalities{m}), 'L')
        BF.inverse.(modalities{m}).L = S.L;
    end
end

bf_save(BF);

out.BF{1} = fullfile(outdir, 'BF.mat');
end


function dep = bf_inverse_vout(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'BF.mat file';
% reference field "B" from output
dep.src_output = substruct('.','BF');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});
end
