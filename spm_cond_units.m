function [y,scalefactor] = spm_cond_units(y,n)
% Scales numeric arrays by a multiple of 10^n to avoid numerical overflow
% FORMAT [y,scalefector] = spm_cond_units(y,n)
%   y - y*scalefactor;
%   n - default 3
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_cond_units.m 4768 2012-06-11 17:06:55Z karl $
 
% default n = 1
%--------------------------------------------------------------------------
try, n; catch, n = 1; end

switch lower(n)
    
    case{'csd'}
        
        % normalise to total power
        %------------------------------------------------------------------
        g     = 0;
        for i = 1:length(y)
            csd = y{i};
            for i = 1:size(csd,2);
                g = g + sum(csd(:,i,i));
            end
        end
        scalefactor = (length(y)*size(csd,1)*size(csd,2))/g;
        y           = spm_unvec(spm_vec(y)*scalefactor,y);
        
    otherwise
        
        % rescale
        %------------------------------------------------------------------
        d           = spm_vec(y);
        scalefactor = std(d(~isnan(d)));
        scalefactor = (10^n)^-round(log10(scalefactor)/n);
        y           = spm_unvec(d*scalefactor,y);
end
