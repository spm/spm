function [y dydp] = spm_gen_erp_dxdp(P,M,U)
% Generates a prediction of trial-specific source activity
% FORMAT [y dydp] = spm_gen_erp_dxdp(P,M,U)
%
% P - parameters
% M - neural-mass model structure
% U - trial-specific effects
%
% y - {[ns,nx],...} - predictions for nx states {trials}
%                   - for ns samples
% dydp              - dertivative w.r.t. p
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_gen_erp_dxdp.m 2353 2008-10-17 11:56:15Z karl $

% within-trial inputs
%==========================================================================
persistent x

% check input u = f(t,P,M)
%--------------------------------------------------------------------------
try, fu  = M.fu; catch,  fu  = 'spm_erp_u'; end
try, ns  = M.ns; catch,  ns  = 128;         end
try, dt  = U.dt; catch,  dt  = 0.004;       end

% peri-stimulus time inputs
%--------------------------------------------------------------------------
U.u = feval(fu,[1:ns]*U.dt,P,M);


% between-trial inputs
%==========================================================================
try
    X = U.X;
    nu = size(U.u,2);
catch
    X = sparse(1,0);
    nu = 0;
end


% cycle over trials
%--------------------------------------------------------------------------
y      = cell(size(X,1),1);
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
    [x dxdp]  = spm_int_LP(Q,M,U);
    y{c}      = x;
    dydp(c,:) = dxdp;

end

% concatenate for spm_nlsi_N
%--------------------------------------------------------------------------
y  = spm_cat(y,1);
df = spm_cell_swap(df);
df = spm_cat(df,1);
