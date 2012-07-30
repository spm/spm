function [f] = spm_fy_erp(y,M)
% feature selection for erp models 
% FORMAT [f] = spm_fy_erp(y,M)
% f = y*M.U;
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_fy_erp.m 4814 2012-07-30 19:56:05Z karl $

% projectors
%--------------------------------------------------------------------------
try, M.U; catch, M.U = 1; end      % spatial
try, M.S; catch, M.S = 1; end      % temporal

% spatial (E) and temporal (S) projection
%--------------------------------------------------------------------------
if isnumeric(y)
    try
        f = M.S'*y*M.U;
    catch
        f = M.U*y*M.S';
    end
else
    for i = 1:length(y)
        f{i} = spm_fy_erp(y{i},M);
    end
end
