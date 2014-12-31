function [f] = spm_fy_erp(y,M)
% feature selection for erp models 
% FORMAT [f] = spm_fy_erp(y,M)
% f = y*M.U;
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_fy_erp.m 6294 2014-12-31 16:47:47Z karl $

% projectors
%--------------------------------------------------------------------------
try, M.U; catch, M.U = 1; end      % spatial

% spatial projection
%--------------------------------------------------------------------------
if isnumeric(y)
      f = y*M.U;
else
    for i = 1:length(y)
        f{i} = spm_fy_erp(y{i},M);
    end
end
