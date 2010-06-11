%% List of modules
% Each module batch is assumed to contain exactly one cfg_exbranch, and
% this should be the last output of the batch
cfg_dir = '/data/projects/spm-devel/matlabbatch/branches/mlb_1/cfg_basicio';
modules = cfg_getfile('fplist', fullfile(cfg_dir,'src'),'^batch_basicio_[0-9]*_.*\.m$');
%% Top level batch
% This batch needs to be modified if the number of modules changes - the
% top level repeat will have to contain the correct number of value
% entries.
toplevel = fullfile(cfg_dir,'src','batch_basicio_toplevel.m');
%% Get modules
exbranches = cell(size(modules));
for cm = 1:numel(modules)
    id = cfg_util('initjob',modules{cm});
    cfg_util('run',id);
    o = cfg_util('getalloutputs',id);
    exbranches{cm} = o{end};
    cfg_util('deljob',id);
end
%% Assemble top level
id = cfg_util('initjob',toplevel);
cfg_util('filljob',id,exbranches{:},{cfg_dir});
cfg_util('run',id);
cfg_util('deljob',id);