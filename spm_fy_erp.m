function [f] = spm_fy_erp(y,M)
% feature selection for erp models 
% FORMAT [f] = spm_fy_erp(y,M)
% f = y*M.E;
%__________________________________________________________________________

% projectors
%--------------------------------------------------------------------------
try, M.E; catch, M.E = 1; end      % spatial
try, M.S; catch, M.S = 1; end      % temporal

% projection
%--------------------------------------------------------------------------
try
    f = M.S'*y*M.E;
catch
    f = M.E*y*M.S';
end