function f = spm_fy_erp(y,M)
% Feature selection for erp models 
% FORMAT f = spm_fy_erp(y,M)
% f = y*M.U;
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging

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
