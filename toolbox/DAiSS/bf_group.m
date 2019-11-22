function out = bf_group
% A module for applying a processing step to a group of subjects
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_group.m 7703 2019-11-22 12:06:29Z guillaume $

%--------------------------------------------------------------------------
% BF
%--------------------------------------------------------------------------
BF = cfg_files;
BF.tag = 'BF';
BF.name = 'BF.mat or M/EEG files';
BF.filter = '.*.mat$';
BF.num = [1 Inf];
BF.help = {'Select BF.mat file.'};


%--------------------------------------------------------------------------
% prefix
%--------------------------------------------------------------------------
prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Directory prefix';
prefix.help    = {'Specify the string to be prepended to the output directory name',...
    'This is only used for batches starting with the Data module.'};
prefix.strtype = 's';
prefix.num     = [0 Inf];
prefix.val     = {''};


%--------------------------------------------------------------------------
% plugin
%--------------------------------------------------------------------------
plugin      = cfg_choice;
plugin.tag  = 'plugin';
plugin.name = 'Group analysis ';

group_funs = spm_select('List', fileparts(mfilename('fullpath')), '^bf_group_.*\.m$');
group_funs = cellstr(group_funs );
for i = 1:numel(group_funs)
    plugin.values{i} = feval(spm_file(group_funs{i},'basename'));
end

out = cfg_exbranch;
out.tag = 'group';
out.name = 'Group analysis';
out.val = {BF, prefix, plugin};
out.help = {'Set up group analyses'};
out.prog = @bf_group_run;
out.vout = @bf_group_vout;
out.modality = {'EEG'};
end

function  out = bf_group_run(job)


plugin_name = cell2mat(fieldnames(job.plugin));

BF          = job.BF;
S           = job.plugin.(plugin_name);
S.prefix    = job.prefix;

BF = feval(['bf_group_' plugin_name], BF, S);

if ~isa(BF, 'cell')
    BF = {BF};
end

out.BF = BF(:);

end

function dep = bf_group_vout(job)
% Output is always in field "BF", no matter how job is structured
dep = cfg_dep;
dep.sname = 'BF.mat files';
% reference field "B" from output
dep.src_output = substruct('.','BF');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end
