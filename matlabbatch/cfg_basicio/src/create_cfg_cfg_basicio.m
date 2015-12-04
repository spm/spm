function create_cfg_cfg_basicio
src_dir = fileparts(mfilename('fullpath'));
cfg_dir = fileparts(src_dir);
out = process_dir(src_dir);
%% Code Generator batch
codegen = fullfile(src_dir,'batch_basicio_codegen.m');
id = cfg_util('initjob',codegen);
cfg_util('filljob',id,{cfg_dir},out);
cfg_util('run',id);
cfg_util('deljob',id);

function out = process_dir(d)
[m, sd] = cfg_getfile('fplist', d,'^batch_basicio_[0-9]*_.*\.m$');
%% List of modules
% Each module batch is assumed to contain exactly one cfg_exbranch, and
% this should be the last output of the batch
outm = cell(size(m));
for cm = 1:numel(m)
    id = cfg_util('initjob',m{cm});
    cfg_util('run',id);
    o = cfg_util('getalloutputs',id);
    outm{cm} = o{end};
    cfg_util('deljob',id);
end
%% List of subdirs
[u, n, e] = cellfun(@fileparts,sd,'UniformOutput',false);
sd = sd(cellfun(@isempty,regexp(strcat(n,e),'^\.')));
outs = cell(size(sd));
for cs = 1:numel(sd)
    outs{cs} = process_dir(sd{cs});
end
outa = [outm(:); outs(:)];
%% Top level batch
% This batch needs to be modified if the number of modules or sublevels
% changes - the top level choice will have to contain the correct number of
% value entries.
if isempty(outa)
    out = {};
else
    fprintf('Processing %s\n', d);
    toplevel = fullfile(d,'batch_basicio_toplevel.m');
    id = cfg_util('initjob',toplevel);
    cfg_util('filljob',id,outa{:});
    cfg_util('run',id);
    o = cfg_util('getalloutputs',id);
    out = o{end};
    cfg_util('deljob',id);
end