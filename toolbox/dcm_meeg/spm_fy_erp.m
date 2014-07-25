function [f] = spm_fy_erp(y,M)
% feature selection for erp models 
% FORMAT [f] = spm_fy_erp(y,M)
% f = y*M.U;
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_fy_erp.m 6123 2014-07-25 17:10:51Z karl $

% projectors
%--------------------------------------------------------------------------
try, M.U; catch, M.U = 1; end      % spatial
try, M.R; catch, M.R = 1; end      % temporal

% spatial (E) and temporal (S) projection
%--------------------------------------------------------------------------
if isnumeric(y)
%     n = size(y,1)/size(M.R,2);
%     f = kron(eye(n,n),M.R)*;
      f = y*M.U;
else
    for i = 1:length(y)
        f{i} = spm_fy_erp(y{i},M);
    end
end
