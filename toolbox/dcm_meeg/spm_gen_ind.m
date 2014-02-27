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
% $Id: spm_gen_ind.m 5900 2014-02-27 21:54:51Z karl $


% within-trial inputs
%==========================================================================

% check input u = f(t,P,M)
%--------------------------------------------------------------------------
try
    fu  = M.fu;
catch
    fu  = 'spm_erp_u';
end

% peri-stimulus time inputs
%--------------------------------------------------------------------------
t   = (1:M.ns)*U.dt;
U.u = feval(fu,t,P,M);


% between-trial inputs
%==========================================================================
try, X = U.X; catch, X = sparse(1,0); end


% cycle over trials
%--------------------------------------------------------------------------
nX     = size(X,1);
y      = cell(nX,1);
for  c = 1:nX
    
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
    y{c}  = spm_int_L(Q,M,U);

end


