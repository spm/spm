function [BF, matlabbatch, view] = bf_wizard_view(S)
% A handy command-line based batch filler with some defaults for DAiSS
% view module, pick a few options, and it will default for unpopulated
% fields
%
% Currently supported output methods include:
%   - glass
%   - surface
%__________________________________________________________________________

% George O'Neill
% Copyright (C) 2022-2023 Wellcome Centre for Human Neuroimaging


if ~isfield(S,'batch'), matlabbatch = []; else; matlabbatch = S.batch;  end
if ~isfield(S,'BF'),        error('I need a BF.mat file specified!');   end
if ~isfield(S,'method'),    error('You need to specify a method!');     end
if ~isfield(S,S.method),    S.(S.method) = struct();                    end
if ~isfield(S,'run'),       S.run = 1;                                  end


% specify BF, ensure its a cell...
if ~iscell(S.BF)
    S.BF = {S.BF};
end

view = struct();
view.BF = S.BF;

try
    opts = feval(['bf_view_' S.method]);
catch
    error('not a valid view method!')
end

view.plugin.(S.method) = struct();

for ii = 1:numel(opts.val)
    
    tag = opts.val{ii}.tag;
    val = opts.val{ii}.val;
    
    if ~isfield(S.(S.method),tag)
        view.plugin.(S.method).(tag) = val{1};
    else
        view.plugin.(S.method).(tag) = S.(S.method).(tag);
    end
     
end

% determine number of jobs in list then iterate by 1;
jobID = numel(matlabbatch) + 1;
% generate matlabbatch
matlabbatch{jobID}.spm.tools.beamforming.view = view;

% Run job (if required)
if S.run
    out = spm_jobman('run',matlabbatch);
    BF = out{1,1}.BF{:};
else
    BF = [];
end