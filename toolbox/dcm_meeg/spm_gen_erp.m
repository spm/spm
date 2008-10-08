function [y,df] = spm_gen_erp(P,M,U)
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
% $Id: spm_gen_erp.m 2315 2008-10-08 14:43:18Z jean $

if ~isfield(M,'fastEM')
    M.fastEM = 0;
end


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
t   = [1:ns]*U.dt;
if ~M.fastEM
    U.u = feval(fu,t,P,M);
else
    [U.u,U.dudP] = feval(fu,t,P,M);
end



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
y = cell(size(X,1),1);

if M.fastEM
    df = cell(size(X,1),1);
    n = size(P.A{1},1);
    In = eye(n);
    indA = spm_getvec(P,'A');
    indH = spm_getvec(P,'H');
    dfdH = sparse(M.n*M.ns,n.^2);
    dfdH(:,In==1) = 1;
    ind2 = find(dfdH);
end

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
        if ~M.fastEM
            y{c}      = spm_int_L(Q,M,U);
            df          = [];
        else
            [y{c},df{c}]   = spm_int_df(Q,M,U);
            % Insert gradients wrt B{j} from gradients wrt A and H...
            for j = 1:size(X,2)
                [indBj] = spm_getvec(Q,'B',j);
                dfdA = df{c}(:,indA);
                tmp = df{c}(:,indH);
                dfdH(ind2) = tmp(:);
                dfdB = X(c,j)*[dfdA,dfdH]*kron(ones(4,1),speye(n*n));
                df{c} = [df{c}(:,1:indBj(1)-1),dfdB,...
                    df{c}(:,indBj(1):end)];
            end
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

