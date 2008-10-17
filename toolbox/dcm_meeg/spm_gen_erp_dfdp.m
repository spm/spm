function [y,df] = spm_gen_erp_dxdp(P,M,U)
% Generates a prediction of trial-specific source activity
% FORMAT [x,dxdp] = spm_gen_erp_dxdp(P,M,U)
%
% P - parameters
% M - neural-mass model structure
% U - trial-specific effects
%
% x - [ns*nc,nx]  - predictions for nx states {trials}
%                 - for ns samples and nc conditions
% dxdp - {1 x np} - derivative of prediction with respect to parameters
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_gen_erp_dfdp.m 2353 2008-10-17 11:56:15Z karl $
 
% within-trial inputs
%==========================================================================
persistent x
try, M.loop; catch, M.loop = 0; end
A   = 1/2;
 
% check input u = f(t,P,M)
%--------------------------------------------------------------------------
try, ns  = M.ns; catch,  ns  = 128;         end
try, dt  = U.dt; catch,  dt  = 0.004;       end
 
% peri-stimulus time inputs
%--------------------------------------------------------------------------
t            = [1:ns]*U.dt;
[U.u,U.dudP] = spm_erp_u_dudp(t,P,M);
 
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
y    = cell(size(X,1),1);
df   = cell(size(X,1),1);
n    = size(P.A{1},1);
In   = eye(n);
indA = spm_getvec(P,'A');
indH = spm_getvec(P,'H');
dfdH = sparse(M.n*M.ns,n.^2);
dfdH(:,In==1) = 1;
ind2 = find(dfdH);
 
 
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
 
        % trial-specific initial states
        %------------------------------------------------------------------
        try, M.x = x{c}; end
 
        % integrate & insert gradients wrt B{j} from gradients wrt A and H...
        %------------------------------------------------------------------
        [y{c},df{c}]   = spm_int_df(Q,M,U);
 
        for j = 1:size(X,2)
            [indBj]    = spm_getvec(Q,'B',j);
            dfdA       = df{c}(:,indA);
            tmp        = df{c}(:,indH);
            dfdH(ind2) = tmp(:);
            dfdB       = X(c,j)*[dfdA,dfdH]*kron(ones(4,1),speye(n*n));
            df{c}      = [df{c}(:,1:indBj(1)-1),dfdB,df{c}(:,indBj(1):end)];
        end
 
 
 
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
 
% rearrange predictions and derivatives in the form expected by spm_nlsi_N
%--------------------------------------------------------------------------
y  = spm_cat(y,1);
df = spm_cell_swap(df);
df = spm_cat(df,1);

