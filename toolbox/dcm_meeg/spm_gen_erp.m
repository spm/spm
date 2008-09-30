function [y] = spm_gen_erp(P,M,U)
% Generates a prediction of trial-specific source activity
% FORMAT [y] = spm_gen_erp(P,M,U)
%
% P - parameters
% M - neural-mass model structure
% U - trial-specific effects
%
% y - {[ns,nx],...} - predictions for nx states {trials}
%                   - for ns samples
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_gen_erp.m 2265 2008-09-30 19:22:48Z karl $


% within-trial inputs
%==========================================================================
persistent x
try, M.loop; catch, M.loop = 0; end
A   = 1/2;

% check input u = f(t,P,M)
%--------------------------------------------------------------------------
try, fu  = M.fu; catch,  fu  = 'spm_erp_u'; end
try, ns  = M.ns; catch,  ns  = 128;         end
try, dt  = U.dt; catch,  dt  = 0.004;       end

% peri-stimulus time inputs
%--------------------------------------------------------------------------
U.u = feval(fu,[1:ns]*dt,P,M);


% between-trial inputs
%==========================================================================
try, X   = U.X; catch,   X   = sparse(1,0); end


% cycle over trials
%--------------------------------------------------------------------------
for  c = 1:size(X,1)

    % baseline parameters
    %----------------------------------------------------------------------
    Q  = P;

    % trial-specific effects
    %----------------------------------------------------------------------
    for i = 1:size(X,2)

        Q.A{1}  = Q.A{1} + X(c,i)*P.B{i};         % forward   connections
        Q.A{2}  = Q.A{2} + X(c,i)*P.B{i};         % backward  connections
        Q.A{3}  = Q.A{3} + X(c,i)*P.B{i};         % lateral   connections

        try
            Q.H = Q.H + X(c,i)*diag(P.B{i});      % intrinsic connections
        catch
            Q.G = Q.G + X(c,i)*diag(P.B{i});
        end
    end


    % integrate DCM for this trial
    %----------------------------------------------------------------------
    for i = 1:4

        % trial-specific intial states
        %------------------------------------------------------------------
        try, M.x = x{c}; end

        % integrate
        %------------------------------------------------------------------
        y{c,1}   = spm_int_L(Q,M,U);

        % save the final state vector
        %------------------------------------------------------------------
        if M.loop
            x{c} = A*spm_vec(y{c}(end,:)) + (1 - A)*spm_vec(M.x);
            x{c} = spm_unvec(x{c},M.x);
        else
            x{c} = M.x; break
        end
    end
end

