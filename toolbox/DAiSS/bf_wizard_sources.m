function [matlabbatch, sources] = bf_wizard_sources(S)

% A handy command-line based batch filler with some defaults for DAiSS
% source module, pick a few options, and it will default for unpopulated
% fields

if ~isfield(S,'batch'); matlabbatch = []; else; matlabbatch = S.batch;  end
if ~isfield(S,'BF');        error('I need a BF.mat file specified!');   end
if ~isfield(S,'reduce_rank'); S.reduce_rank = [2 3];                    end
if ~isfield(S,'keep3d');    S.keep3d = true;                            end
if ~isfield(S,'method');    error('You need to specify a method!');     end
if ~isfield(S,S.method);    S.(S.method) = struct();                    end
if ~isfield(S,'normalise_lf'); S.normalise_lf = false;                  end
if ~isfield(S,'visualise'); S.visualise = true;                         end


% specify BF, ensure its a cell...
if ~iscell(S.BF)
    S.BF = {S.BF};
end

sources = struct();
sources.BF = S.BF;
sources.reduce_rank = S.reduce_rank;


% determine if the source space method is valid
try
    opts = feval(['bf_sources_' S.method]);
    % check its only either the mesh or grid method for now
    % TODO: test wizard support for phantom grid, mni coord or VOIs
    target = {'mesh','grid','scalp','mni_coords'};
    if ~ismember(S.method,target)
        error(['method ' S.method ' is yet to be supported through the '...
            'wizard, please use another [mesh, grid, scalp]']);
    end
catch
    error('not a valid source space generation method!')
end

% Special cases
switch S.method
    case 'mni_coords'
        if ~isfield(S.(S.method),'pos')
            error('please specify some mni coordinates!')
        end
end

% Poplate some plugin options
sources.plugin.(S.method) = struct();

for ii = 1:numel(opts.val)
    
    tag = opts.val{ii}.tag;
    val = opts.val{ii}.val;
    
    if ~isfield(S.(S.method),tag)
        sources.plugin.(S.method).(tag) = val{1};
    else
        sources.plugin.(S.method).(tag) = S.(S.method).(tag);
    end
     
end

% Final options
sources.normalise_lf = S.normalise_lf;
sources.visualise = S.visualise;

% determine number of jobs in list then iterate by 1;
jobID = numel(matlabbatch) + 1;
% generate matlabbatch
matlabbatch{jobID}.spm.tools.beamforming.sources = sources;

