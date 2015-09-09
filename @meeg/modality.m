function [res, list] = modality(this, scalp, planar)
% Returns data modality 
% FORMAT [res, list] = modality(this, scalp)
%
% scalp - 1 (default) only look at scalp modalities
%         0  look at all modalities
% planar - 1 distinguish between MEG planar and other MEG 
%          0 (default) do not distinguish
% If more than one modality is found the function returns 'Multimodal'
% in res and a cell array of modalities in list.
% _______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: modality.m 6542 2015-09-09 11:48:34Z karl $

if nargin == 1
    scalp = 1;
end
if nargin < 3
    planar = 0;
end

list = {};

if ~isempty(indchantype(this, {'MEG', 'MEGPLANAR', 'MEGCOMB'}))
    if planar
        if ~isempty(indchantype(this, 'MEGPLANAR'))
            list = [list {'MEGPLANAR'}];
        end
        if ~isempty(indchantype(this, 'MEG'))
            list = [list {'MEG'}];
        end
        if ~isempty(indchantype(this, 'MEGCOMB'))
            list = [list {'MEGCOMB'}];
        end
    else
        list = [list {'MEG'}];
    end
end

if ~isempty(indchantype(this, 'EEG'))
    list = [list {'EEG'}];
end

if  ~isempty(indchantype(this, 'LFP')) && ~scalp
    list = [list {'LFP'}];
end

if  ~isempty(indchantype(this, 'ILAM')) && ~scalp
    list = [list {'ILAM'}];
end

switch numel(list)
    case 0
        res = 'Other';
    case 1
        res = list{1};
    otherwise
        res = 'Multimodal';
end