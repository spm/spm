function [y,scale] = spm_cond_units(y,n)
% Scales numeric arrays by a multiple of 10^n to avoid numerical overflow
% FORMAT [y,scale] = spm_cond_units(y,n)
%   y - y*scale;
%   n - default 3
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_cond_units.m 6110 2014-07-21 09:36:13Z karl $
 
% default n = 1
%--------------------------------------------------------------------------
if nargin < 2, n = 1; end

switch lower(n)
    
    case{'csd'}
        
        % normalise to total power
        %------------------------------------------------------------------
        g     = 0;
        for i = 1:length(y)
            csd = y{i};
            for j = 1:size(csd,2);
                g = g + sum(csd(:,j,j));
            end
        end
        scale = (length(y)*size(csd,1)*size(csd,2))/g;
        y     = spm_unvec(spm_vec(y)*scale,y);
        
    otherwise
        
        % rescale
        %------------------------------------------------------------------
        if isnumeric(y)
            scale = max(max(max(abs(y))))/8;
            scale = (10^n)^-round(log10(scale)/n);
            y     = y*scale;
        else
            d     = spm_vec(y);
            scale = max(abs(d))/8;
            scale = (10^n)^-round(log10(scale)/n);
            y     = spm_unvec(d*scale,y);
        end
end
