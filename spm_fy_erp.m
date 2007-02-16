function [f] = spm_fy_erp(y,M)
% feature selection for erp models 
% FORMAT [f] = spm_fy_erp(y,M)
% f = y*M.E;
%__________________________________________________________________________

% projection
%--------------------------------------------------------------------------
try
    f = y*M.E;
catch
    f = M.E'*y;
end
