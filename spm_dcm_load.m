function DCM = spm_dcm_load(P)
% Loads a cell array of DCM filenames into a subjects x models cell array
% FORMAT DCM = spm_dcm_load(P)
%
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% Peter Zeidman
% $Id$

if ~iscell(P)
    P = {P};
end

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
                error('Failed to load model for subject %s model %m', s, m);
            end                
        end
    end
end