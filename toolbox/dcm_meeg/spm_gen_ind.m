function [y] = spm_gen_ind(P,M,U)
% Generates a prediction of trial-specific induced activity
% FORMAT [y] = spm_gen_ind(P,M,U)
%
% P - parameters
% M - neural-mass model structure
% U - trial-specific effects
%
% y - prediction
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_gen_ind.m 1208 2008-03-13 20:59:12Z karl $


% within-trial inputs
%==========================================================================

% check input u = f(t,P,M)
%--------------------------------------------------------------------------
try
    fu  = M.fu;
catch
    fu  = 'spm_ind_u';
end

% peri-stimulus time inputs
%--------------------------------------------------------------------------
t   = [1:M.ns]*U.dt;
U.u = feval(fu,t,P,M);


% between-trial inputs
%==========================================================================
try, X = U.X; catch, X = sparse(1,0); end


% cycle over trials
%--------------------------------------------------------------------------
for  c = 1:size(X,1)
    
    % baseline parameters
    %----------------------------------------------------------------------
    Q  = P;

    % trial-specific inputs
    %----------------------------------------------------------------------
    for i = 1:size(X,2)
          Q.A = Q.A + X(c,i)*Q.B{i};
    end

    % integrate DCM for this trial
    %----------------------------------------------------------------------
    y{c,1} = spm_int_L(Q,M,U);

end


