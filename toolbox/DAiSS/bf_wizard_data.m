function [matlabbatch, data] = bf_wizard_data(S)

% A handy command-line based batch filler with some defaults for DAiSS
% data module, pick a few options, and it will default for unpopulated
% fields
% FORMAT [batch, data] = bf_wizard_data(S)
%   S               - input structure
% Optional fields of S:
%   S.D             - SPM MEEG object               - Default: REQUIRED
%   S.dir           - path to save DAiSS BF.mat     - Default: same as S.D
%   S.val           - which D.inv to use            - Default: 1
%   S.gradsource    - where to pool sensor information from 
%                       (inv | sens)                
%                                                   - Default: 'inv'
%   S.space         - which space to do calculations in
%                       (MNI-Aligned | Head | Native)
%                                                   - Default: MNI-Aligned
%   S.overwite      - Overwrite existing BF.mat     - Default: 0

% Output:
%  matlabbatch      - matlabbatch job for spm_jobman to run
%  data             - simplified summary of options selected
%__________________________________________________________________________
% Copyright (C) 2022 Wellcome Centre for Human Neuroimaging

if ~isfield(S,'batch'); matlabbatch = []; else; matlabbatch = S.batch;  end
if ~isfield(S,'dir');   S.dir = [];                                     end
if ~isfield(S,'D');     error('I need a SPM MEEG dataset provided!');   end
if ~isfield(S,'val');   S.val = 1;                                      end
if ~isfield(S,'gradsource'); S.gradsource = 'inv';                      end
if ~isfield(S,'space'); S.space = 'MNI-aligned';                        end
if ~isfield(S,'overwrite'); S.overwrite = 0;                            end

% Check if dir is empty and determine path of S.D
if isempty(S.dir)
    if isa(S.D,'meeg')
        S.dir = path(S.D);
    else
        D = spm_eeg_load(S.D);
        S.dir =  path(S.D);
    end
    warning(['No explicit path defined, '...
        'placing in same location as dataset']);
end

% Check if the gradsource is correctly defined
target = {'inv','sens'};
if ~ismember(S.gradsource,target)
    error('Please specify a valid gradsource [inv, sens]');
end
    
% Check if the space is correctly defined
target = {'MNI-aligned', 'Head', 'Native'};
if ~ismember(S.space,target)
    error('Please specify a valid space [MNI-aligned, Head, Native]');
end

% determine number of jobs in list then iterate by 1;
jobID = numel(matlabbatch) + 1;
% generate matlabbatch
data = S;
if isfield(data,'batch')
    data = rmfield(data,'batch');
end
if ~iscell(data.dir)
    data.dir = {data.dir};
end
if isa(data.D,'meeg')
    data.D = fullfile(data.D.path,data.D.fname);
end
if ~iscell(data.D)
    data.D = {data.D};
end
matlabbatch{jobID}.spm.tools.beamforming.data = data;
