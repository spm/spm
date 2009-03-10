function D = spm_eeg_sort_conditions(S)
% Define custom order for conditions in an M/EEG file
% FORMAT D = spm_eeg_sort_conditions(S)
% S - existing configuration struct (optional)
% Fields of S:
%   S.D - file name or SPM MEEG dataset
%   S.condlist - cell array of strings which is a subset of condition
%                labels in the file.
%   S.save - 1 - save the dataset (default)
%            0 - do not save
% OUTPUT:
%   D - modified dataset
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_sort_conditions.m 2850 2009-03-10 21:54:38Z guillaume $

SVNrev = '$Rev: 2850 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG Condition sorting'); spm('Pointer','Watch');

%-Get MEEG object
%--------------------------------------------------------------------------
try
    D = S.D;
catch
    [D, sts] = spm_select(1, 'mat', 'Select M/EEG mat file');
    if ~sts, D = []; return; end
    S.D = D;
end

D = spm_eeg_load(D);

%-Do the sorting
%--------------------------------------------------------------------------
if ~isfield(S, 'condlist')
    oldcondlist = D.condlist;
    S.condlist = cell(size(oldcondlist));
    for i = 1:D.nconditions
        str = sprintf('%s|', oldcondlist{:});
        str = str(1:(end-1));

        ind = spm_input(['Select condition ' num2str(i)], 1, 'm', str, 1:numel(oldcondlist));
        S.condlist(i) = oldcondlist(ind); 
        oldcondlist(ind) = [];
    end
end

D = condlist(D, S.condlist);

D = history(D, 'spm_eeg_sort_conditions', S);

if ~isfield(S, 'save') || S.save
    save(D);
end

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG Condition sorting: done'); spm('Pointer','Arrow');
