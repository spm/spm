function [matlabbatch, write] = bf_wizard_write(S)

if ~isfield(S,'batch'); matlabbatch = []; else; matlabbatch = S.batch; end
if ~isfield(S,'BF'); error('I need a BF.mat file specified!'); end
if ~isfield(S,'method'); error('You need to specify an output method!'); end

% specify BF, ensure its a cell...
if ~iscell(S.BF)
    S.BF = {S.BF};
end

write = struct();
write.BF = S.BF;

% catch some common terms for outputs and name them correctly
switch lower(S.method)
    case {'vol','volume','volumetric'}
        S.method_orig  = S.method;
        S.method = 'nifti';
    case {'surf','surface','mesh'}
        S.method_orig  = S.method;
        S.method = 'gifti';
    case {'d','spm','meeg','dataset','spmmeeg'}
        S.method_orig = S.method;
        S.method = 'spmeeg';
end

if isfield(S,'method_orig') % incase a rename has happened;
    if isfield(S,S.method_orig)
        S.(S.method) = S.(S.method_orig);
    end
end

% final catch for later
if ~isfield(S,S.method); S.(S.method) = struct(); end

% special defaults which cant be automatically parsed
switch S.method
    case 'spmeeg'
        if ~isfield(S.spmeeg,'addchannels')
            S.spmeeg.addchannels = 'none';
        end
end

% pull some options from cfg
try
    opts = feval(['bf_write_' S.method]);
catch
    error('not a valid write method!')
end

write.plugin.(S.method) = struct();

for ii = 1:numel(opts.val)
    
    tag = opts.val{ii}.tag;
    val = opts.val{ii}.val;
    
    if ~isfield(S.(S.method),tag)
        write.plugin.(S.method).(tag) = val{1};
    else
        write.plugin.(S.method).(tag) = S.(S.method).(tag);
    end
    
end

% determine number of jobs in list then iterate by 1;
jobID = numel(matlabbatch) + 1;
% generate matlabbatch
matlabbatch{jobID}.spm.tools.beamforming.write = write;