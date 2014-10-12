function [csd] = spm_csd_sel(csd,M)
% Time frequency response - feature selection
% FORMAT [csd] = spm_csd_sel(csd,M)
%
% csd  - {Y(t,w,nc,nc}} - cross-spectral density for nc channels {trials}
%                       - for w frequencies over time t in M.Hz
% M - neural mass model structure
%__________________________________________________________________________
%
% Feature selection for time-ftequency DCMs
%__________________________________________________________________________
% Copyright (C) 2012-2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_csd_sel.m 6234 2014-10-12 09:59:10Z karl $


% baseline correct
%--------------------------------------------------------------------------
for c = 1:numel(csd)
    for w = 1:size(csd{c},2)
        for i = 1:size(csd{c},3)
            for j = 1:size(csd{c},4)
                csd{c}(:,w,i,j) = csd{c}(:,w,i,j) - mean(real(csd{c}(:,w,i,j)));
            end
        end
    end  
end