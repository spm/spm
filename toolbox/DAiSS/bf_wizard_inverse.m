function [BF, matlabbatch, inverse] = bf_wizard_inverse(S)
% A handy command-line based batch filler with some defaults for DAiSS
% invert module, pick a few options, and it will default for unpopulated
% fields
%__________________________________________________________________________

% George O'Neill
% Copyright (C) 2022-2023 Wellcome Centre for Human Neuroimaging


if ~isfield(S,'batch'), matlabbatch = []; else; matlabbatch = S.batch;  end
if ~isfield(S,'BF'), error('I need a BF.mat file specified!');          end
if ~isfield(S,'method'), error('You need to specify a method!');        end
if ~isfield(S,S.method), S.(S.method) = struct();                       end
if ~isfield(S,'run'), S.run = 1;                                        end


% specify BF, ensure its a cell...
if ~iscell(S.BF)
    S.BF = {S.BF};
end

inverse = struct();
inverse.BF = S.BF;

try
    opts = feval(['bf_inverse_' S.method]);
catch
    error('not a valid inverse method!')
end

inverse.plugin.(S.method) = struct();

for ii = 1:numel(opts.val)

    tag = opts.val{ii}.tag;
    val = opts.val{ii}.val;

    if ~isfield(S.(S.method),tag)
        inverse.plugin.(S.method).(tag) = val{1};
    else
        inverse.plugin.(S.method).(tag) = S.(S.method).(tag);
    end

end

% determine number of jobs in list then iterate by 1;
jobID = numel(matlabbatch) + 1;
% generate matlabbatch
matlabbatch{jobID}.spm.tools.beamforming.inverse = inverse;

if S.run
    out = spm_jobman('run',matlabbatch);
    BF = out{1,1}.BF{:};
else
    BF = [];
end
