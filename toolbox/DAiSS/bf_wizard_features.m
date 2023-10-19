function [BF, matlabbatch, features] = bf_wizard_features(S)
% A handy command-line based batch filler with some defaults for DAiSS
% features module, pick a few options, and it will default for unpopulated
% fields
%
% FORMAT [BF, batch, features] = bf_wizard_data(S)
%   S               - input structure
%
% Output:
%  BF               - Resultant DAiSS BF structure
%  batch            - matlabbatch job for spm_jobman to run
%  features         - simplified summary of options selected
%__________________________________________________________________________

% George O'Neill
% Copyright (C) 2022-2023 Wellcome Centre for Human Neuroimaging


if ~isfield(S,'batch'), matlabbatch = []; else; matlabbatch = S.batch;  end
if ~isfield(S,'BF'), error('I need a BF.mat file specified!');          end
if ~isfield(S,'conditions'),    S.conditions = 'all';                   end
if ~isfield(S,'modality'),      S.modality = {'MEG'};                   end
if ~isfield(S,'fuse'),          S.fuse = 'no';                          end
if ~isfield(S,'cross_terms'),   S.cross_terms = 'megeeg';               end
if ~isfield(S,'woi'),           S.woi = [-Inf Inf];                     end
if ~isfield(S,'method'),        S.method = 'identity';                  end
if ~isfield(S,'reg'),           S.reg = 'none';                         end
if ~isfield(S,'bootstrap'),     S.bootstrap = false;                    end
if ~isfield(S,'visualise'),     S.visualise = false;                    end
if ~isfield(S,S.method),        S.(S.method) = struct();                end
if ~isfield(S,S.reg),           S.(S.reg) = struct();                   end
if ~isfield(S,'run'),           S.run = 1;                              end

% specify BF, ensure its a cell...
if ~iscell(S.BF)
    S.BF = {S.BF};
end

features = struct();
features.BF = S.BF;

% check the response to conditions
if sum(strcmp(S.conditions,'all'))
    features.whatconditions.all = 1;
else
    % check its a cell aray
    if ~iscell(S.conditions)
        error('custom condition labels must be in a cell array!')
    end
    features.whatconditions.condlabel = S.conditions;
end

if ~iscell(S.modality)
    S.modality = {S.modality};
end

features.modality = S.modality;
features.fuse = S.fuse;
features.cross_terms = S.cross_terms;
features.woi = S.woi;

% Covariance generation section
try
    opts = feval(['bf_features_' S.method]);
catch
    error('not a valid covariance generation method!')
end



features.plugin.(S.method) = struct();

for ii = 1:numel(opts.val)

    tag = opts.val{ii}.tag;
    val = opts.val{ii}.val;

    if ~isfield(S.(S.method),tag)
        features.plugin.(S.method).(tag) = val{1};
    else
        features.plugin.(S.method).(tag) = S.(S.method).(tag);
    end

end

% Regularisation section
if strcmp(S.reg,'none')
    features.regularisation.manual.lambda = 0;
else

    try
        opts = feval(['bf_regularise_' S.reg]);
    catch
        error('not a valid regularisation method!')
    end

    for ii = 1:numel(opts.val)

        tag = opts.val{ii}.tag;
        val = opts.val{ii}.val;

        if ~isfield(S.(S.reg),tag)
            features.regularisation.(S.reg).(tag) = val{1};
        else
            features.regularisation.(S.reg).(tag) = S.(S.reg).(tag);
        end

    end

end

features.bootstrap = S.bootstrap;
features.visualise = S.visualise;

% determine number of jobs in list then iterate by 1;
jobID = numel(matlabbatch) + 1;
% generate matlabbatch
matlabbatch{jobID}.spm.tools.beamforming.features = features;

% run job if required
if S.run
    out = spm_jobman('run',matlabbatch);
    BF = out{1,1}.BF{:};
else
    BF = [];
end
