function batch = spm_reorder_contrasts(SPM,order)
% Recompute contrasts allowing for permutation and deletion
% FORMAT batch = spm_reorder_contrasts(SPM,order)
% SPM   - SPM data structure
% order - array of contrast indices
%
% batch - batch job
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_reorder_contrasts.m 5827 2014-01-02 16:06:52Z guillaume $


%-Get SPM
%--------------------------------------------------------------------------
if ~nargin
    [SPM, sts] = spm_select(1,'^SPM\.mat$','Select SPM.mat');
    if ~sts, batch = []; return; end
end
if ischar(SPM)
    load(SPM);
end
if ~isfield(SPM,'xCon') || isempty(SPM.xCon)
    batch = []; return;
end

%-Get order
%--------------------------------------------------------------------------
if nargin < 2
    order = 1:numel(SPM.xCon);
end

if islogical(order)
    order = find(order);
end

if iscellstr(order)
    [lia,order] = ismember(order,{SPM.xCon.name});
    if ~all(lia)
        warning('Some contrast names were not found.');
        order = order(lia);
    end
end

if any(order > numel(SPM.xCon) | order < 1)
    error('Invalid contrast index.');
end

%-Create contrasts batch job
%--------------------------------------------------------------------------
batch{1}.spm.stats.con.spmmat = cellstr(fullfile(SPM.swd,'SPM.mat'));
batch{1}.spm.stats.con.delete = 1;
for i=1:numel(order)
    switch SPM.xCon(order(i)).STAT
        case 'T'
            batch{1}.spm.stats.con.consess{i}.tcon.name = SPM.xCon(order(i)).name;
            batch{1}.spm.stats.con.consess{i}.tcon.weights = SPM.xCon(order(i)).c';
        case 'F'
            batch{1}.spm.stats.con.consess{i}.fcon.name = SPM.xCon(order(i)).name;
            batch{1}.spm.stats.con.consess{i}.fcon.weights = SPM.xCon(order(i)).c';
        otherwise
            error('Sorry, T&F tests only.');
    end
end

%-Execute or display batch job
%--------------------------------------------------------------------------
if nargin < 2
    spm_jobman('interactive',batch);
else
    spm_jobman('run',batch);
end
