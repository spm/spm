function out = bf_inverse
% Compute inverse projectors
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2015-2023 Wellcome Centre for Human Neuroimaging


out          = cfg_exbranch;
out.tag      = 'inverse';
out.name     = 'Inverse solution';
out.val      = @bf_inverse_cfg;
out.help     = {'Compute inverse projectors'};
out.prog     = @bf_inverse_run;
out.vout     = @bf_inverse_vout;
out.modality = {'EEG'};


%==========================================================================
function varargout = bf_inverse_cfg

persistent cfg
if ~isempty(cfg), varargout = {cfg}; return; end

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

[cfg,varargout{1}] = deal({BF, plugin});


%==========================================================================
function  out = bf_inverse_run(job)

outdir = spm_file(job.BF{1}, 'fpath');

% cd(outdir);

BF = bf_load(fullfile(outdir, 'BF.mat'),...
    {'data', 'sources', 'features'});
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

% bf_save(BF);

bf_save_path(BF,fullfile(outdir, 'BF.mat'));

out.BF{1} = fullfile(outdir, 'BF.mat');


%==========================================================================
function dep = bf_inverse_vout(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'BF.mat file';
% reference field "B" from output
dep.src_output = substruct('.','BF');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});
