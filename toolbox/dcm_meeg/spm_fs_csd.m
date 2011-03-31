function [y] = spm_fs_csd(y,M)
% spectral feature selection for a CSD DCM
% FORMAT [y] = spm_fs_csd(y,M)
% y      - CSD
% M      - model structure
%__________________________________________________________________________
%
% This simply log-transforms the (real) auto-spectra
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fs_csd.m 4281 2011-03-31 19:49:57Z karl $


% just return auto-spectra at present
%--------------------------------------------------------------------------
y     = y;

return

% log-transform (real) auto-spectra
%--------------------------------------------------------------------------
for i = 1:length(y);
    for j = 1:size(y{i},2)
        y{i}(:,j,j) = log(y{i}(:,j,j));
    end
end

