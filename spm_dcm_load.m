function DCM = spm_dcm_load(P)
% Load a cell array of DCM filenames into a subjects x models cell array
% FORMAT DCM = spm_dcm_load(P)
%
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% Peter Zeidman
% $Id: spm_dcm_load.m 6702 2016-01-28 15:10:49Z guillaume $


P   = cellstr(P);
DCM = cell(size(P));

for s = 1:size(P,1)
    for m = 1:size(P,2)
        if isstruct(P{s,m})
            DCM{s,m} = P{s,m};
        else
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
