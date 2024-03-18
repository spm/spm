function [BF,matlabbatch, sources] = bf_wizard_sources(S)
% A handy command-line based batch filler with some defaults for DAiSS
% source module, pick a few options, and it will default for unpopulated
% fields
%
% FORMAT [BF, batch, sources] = bf_wizard_sources(S)
%   S               - input structure
% Optional fields of S:
%   S.BF            - Path to a BF structure or the
%                       structure itself
%                                                   - Default: REQUIRED
%
%   S.batch         - matlabbatch, of which this job
%                       can be appended to
%                                                   - Default: []
%
%   S.reduce_rank   - [1x2] vector determining the 
%                     dimensionality of the lead
%                     fields for a source, first
%                     element for MEG, second for EEG
%                                                   - Default: [2 3]
%
%   S.keep3d        - If the leadfield rank has been
%                     reduced with S.reduce_rank, this
%                     ensures the are still 3 lead fields
%                     per source.               
%                                                   - Default: true
%
%   S.normalise_lf  - Make the norms of each lead field 1
%                                                   - Default: false
% 
%   S.visualise     - Visualise the source space, sensors
%                     and conductive boundar[y/ies]
%                                                   - Default: true
%   
%   S.method        - How do we want the source space
%                     generated? Validated methods with
%                     this are ( 'grid' | 'mesh' )   
%                                                   - Default: REQUIRED
%
%   S.run           - Run the batch, set to 0 to
%                     bypass the run for debugging
%                                                   - Default: 1
%
%
% Method options for S:
% Extensive details can be found at 
% https://github.com/spm/spm/blob/main/toolbox/DAiSS/doc/commands/02_sources.md
% But a summary of the essential options below.
%
% GRID METHOD
%   
%   S.grid.resolution - Distance between sources
%                        (in mm)                    - Default: 5
%
%   S.grid.constrain  - Which boundary do we want 
%                       sources outside of to be
%                       excluded? Options: ('iskull'
%                       | 'oskull' | 'scalp')       - Default: 'iskull'
%
% MESH METHOD
%
%  S.mesh.orient      - How are sources oriented on
%                       the vertices of the mesh?
%                       'unoriented' keeps the lead field
%                       triplet, whilst 'original' returns 
%                       one lead field normal to the mesh
%                       surface
%                                                   - Default: 'unoriented'
% Output:
%  BF               - Resultant DAiSS BF structure
%  batch            - matlabbatch job for spm_jobman to run
%  sources          - simplified summary of options selected
%__________________________________________________________________________

% George O'Neill
% Copyright (C) 2022-2023 Wellcome Centre for Human Neuroimaging


if ~isfield(S,'batch'), matlabbatch = []; else; matlabbatch = S.batch;  end
if ~isfield(S,'BF'),        error('I need a BF.mat file specified!');   end
if ~isfield(S,'reduce_rank'), S.reduce_rank = [2 3];                    end
if ~isfield(S,'keep3d'),    S.keep3d = true;                            end
if ~isfield(S,'method'),    error('You need to specify a method!');     end
if ~isfield(S,S.method),    S.(S.method) = struct();                    end
if ~isfield(S,'normalise_lf'), S.normalise_lf = false;                  end
if ~isfield(S,'visualise'), S.visualise = true;                         end
if ~isfield(S,'run'), S.run = 1;                                        end

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

% Run job (if required)
if S.run
    out = spm_jobman('run',matlabbatch);
    BF = out{1,1}.BF{:};
else
    BF = [];
end


