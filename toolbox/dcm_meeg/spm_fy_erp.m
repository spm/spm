function [f] = spm_fy_erp(y,M)
% feature selection for erp models 
% FORMAT [f] = spm_fy_erp(y,M)
% f = y*M.E;
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_fy_erp.m 1176 2008-02-28 13:29:29Z karl $

% projectors
%--------------------------------------------------------------------------
try, M.E; catch, M.E = 1; end      % spatial
try, M.S; catch, M.S = 1; end      % temporal

% spatial (E) and temporal (S) projection
%--------------------------------------------------------------------------
if isnumeric(y)
    try
        f = M.S'*y*M.E;
    catch
        f = M.E*y*M.S';
    end
else
    for i = 1:length(y)
        try
            f{i} = M.S'*y{i}*M.E;
        catch
            f{i} = M.E*y{i}*M.S';
        end
    end
end
