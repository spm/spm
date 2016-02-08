function DCM = spm_dcm_load(P)
% Load a cell array of DCM filenames into a subjects x models cell array
% FORMAT DCM = spm_dcm_load(P)
%
% P - a character array or cell array of DCM filenames
%
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% Peter Zeidman
% $Id: spm_dcm_load.m 6716 2016-02-08 18:21:37Z peter $


if ischar(P)
    P   = cellstr(P);
end
DCM = cell(size(P));

for s = 1:size(P,1)
    for m = 1:size(P,2)
        if isstruct(P{s,m})
            % DCM already a structure - leave unchanged
            DCM{s,m} = P{s,m};
        else
            % Load DCM from filename
            try
                model    = load(P{s,m});
                DCM{s,m} = model.DCM;
            catch
                fprintf('File: %s\n',P{s,m});
                error('Failed to load model for subject %d model %d', s, m);
            end                
        end
    end
end
