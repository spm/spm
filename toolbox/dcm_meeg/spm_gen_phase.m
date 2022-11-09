function [x] = spm_gen_phase(P,M,U)
% Generate state activities for trial-specific phase-coupled activity
% FORMAT [x] = spm_gen_phase(P,M,U)
%
% P - parameters
% M - model structure
% U - trial-specific effects
%
% x - states
%__________________________________________________________________________

% Will Penny
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


% between-trial inputs
%==========================================================================
try, X = U.X; catch, X = sparse(1,0); end

if isfield(M,'trial')
    Mx=M.x;
end

% cycle over trials
%--------------------------------------------------------------------------
for  c = 1:size(X,1)
    
    % baseline parameters
    %----------------------------------------------------------------------
    Q  = P;

    % trial-specific inputs
    %----------------------------------------------------------------------
    for i = 1:size(X,2)
          Q.As = Q.As + X(c,i)*Q.Bs{i};
    end
    if isfield(Q,'Ac')
        for i = 1:size(X,2)
            Q.Ac = Q.Ac + X(c,i)*Q.Bc{i};
        end
    end
    
    if isfield(M,'trial')
        % Set up trial-specific initial state
        M.x=M.trial{c}.x;
    end
    
    x{c,1} = spm_int_ode(Q,M,U);
    
end

if isfield(M,'trial')
    M.x=Mx;
end

