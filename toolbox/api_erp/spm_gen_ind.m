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
% $Id: spm_gen_ind.m 1174 2008-02-27 20:22:30Z karl $


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

% baseline parameters
%--------------------------------------------------------------------------
Q  = P;

% cycle over trials
%--------------------------------------------------------------------------
for  c = 1:size(X,1)

    % trial-specific inputs
    %----------------------------------------------------------------------
    for i = 1:size(X,2)
          P.A = Q.A + X(c,i)*Q.B{i};
    end

    % integrate DCM for this trial
    %----------------------------------------------------------------------
    y{c,1} = spm_int_L(P,M,U);

end


