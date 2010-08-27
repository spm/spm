function [y,scalefactor] = spm_cond_units(y,n)
% Scales numeric arrays by a multiple of 10^n to avoid numerical overflow
% FORMAT [y,scalefector] = spm_cond_units(y,n)
%   y - y*scalefactor;
%   n - default 3
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_cond_units.m 4052 2010-08-27 19:22:44Z karl $
 
% default n = 3
%--------------------------------------------------------------------------
try, n; catch, n = 3; end

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
        scalefactor = norm(d(~isnan(d)),1);
        scalefactor = (10^n)^-round(log10(scalefactor)/n);
        y           = spm_unvec(d*scalefactor,y);
end
