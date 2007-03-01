function [Ep,Eg,Cp,Cg,S,F] = spm_nlsi_N(M,U,Y)
% Bayesian inversion of a nonlinear model of the form G(g)*F(u,p)
% FORMAT [Ep,Eg,Cp,Cg,S,F]= spm_nlsi_N(M,U,Y)
%
% Dynamical MIMO models
%__________________________________________________________________________
% M.IS - intgration scheme for hidden states f(p,M,U) - generative model
%
%     M.f  - f(x,u,p,M) - state equation:  dxdt = f(x,u)
%     M.G  - G(g,M)     - linear observer: y    = G*x
%
% M.FS - function name f(y,M)   - feature selection
%        This [optional] function perfoms feature selection assuming the
%        generalised model y = FS(y,M) = FS(G*x,M) + X0*P0 + e
%
%
% M.P  - starting estimtes for model parameters [optional]
%
% M.pE - prior expectation  - of model parameters - f(x,u,p,M)
% M.pC - prior covariance   - of model parameters - f(x,u,p,M)
%
% M.gE - prior expectation  - of model parameters - G(g,M)
% M.gC - prior covariance   - of model parameters - G(g,M)
%
% M.hE - prior expectation  - E{h}   of precision parameters
% M.hC - prior covariance   - Cov{h} of precision parameters
%
% U.u  - inputs
% U.dt - sampling interval
%
% Y.y  - outputs
% Y.X0 - Confounds or null space
% Y.dt - sampling interval for outputs
% Y.Q  - error precision components
%
%
% Parameter estimates
%--------------------------------------------------------------------------
% Ep  - (p x 1)         conditional expectation  E{p|y}
% Cp  - (p x p)         conditional covariance   Cov{p|y}
%
% Eg  - (p x 1)         conditional expectation  E{g|y}
% Cg  - (p x p)         conditional covariance   Cov{g|y}
%
% S   - (v x v)         [Re]ML estimate of error Cov{e(h)}
%
% log evidence
%--------------------------------------------------------------------------
% F   - [-ve] free energy F = log evidence = p(y|m)
%
%__________________________________________________________________________
% Returns the moments of the posterior p.d.f. of the parameters of a
% nonlinear model specified by IS(P,M,U) under Gaussian assumptions. Usually,
% IS would be an integrator of a dynamic MIMO input-state-output model 
%
%              dx/dt = f(x,u,p)
%              y     = G(g)*x  + X0*B + e
%
% The E-Step uses a Fisher-Scoring scheme and a Laplace
% approximation to estimate the conditional expectation and covariance of P
% If the free-energy starts to increase, a Levenburg-Marquardt scheme is
% invoked.  The M-Step estimates the precision components of e, in terms
% of [Re]ML point estimators of the log-precisions.
% An optional feature selection can be specified with parameters M.FS
%
%--------------------------------------------------------------------------
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience
 
% Karl Friston
% $Id: spm_nlsi_GN.m 731 2007-02-07 14:31:41Z karl $
 
% figure (unless disabled)
%--------------------------------------------------------------------------
try
    M.nograph;
catch 
    Fsi = spm_figure('GetWin','SI');
    if isempty(Fsi)
       Fsi =  spm_figure('Create','SI','System Identification');
    else
       clf
    end
end

% check integrator
%--------------------------------------------------------------------------
try
    IS = M.IS;
catch
    IS = 'spm_int_U';
end

% check observer has not been accidently specified
%--------------------------------------------------------------------------
try
    M = rmfield(M,'g');
end

% composition of feature selection and prediction (usually an integrator)
%--------------------------------------------------------------------------
if isfield(M,'FS')

    % FS(y,M)
    %----------------------------------------------------------------------
    try
        y  = feval(M.FS,Y.y,M);
        GM = inline([M.FS '(' M.G '(P,M),M)'],'P','M');

    % FS(y)
    %----------------------------------------------------------------------
    catch
        y  = feval(M.FS,Y.y);
        GM = inline([M.FS '(' M.G '(P,M))'],'P','M');
    end
else

    % y
    %----------------------------------------------------------------------
    y   = Y.y;
    GM  = inline([M.G '(P,M)'],'P','M');
end

% data y
%--------------------------------------------------------------------------
[ns nr] = size(y);          % number of samples and responses
M.ns    = ns;               % store in M.ns for integrator

% initial states
%--------------------------------------------------------------------------
try
    M.x;
catch
    if ~isfield(M,'n'), M.n = 0;    end
    M.x = sparse(M.n,1);
end

% input
%--------------------------------------------------------------------------
try
    U;
catch
    U = [];
end

% initial parameters
%--------------------------------------------------------------------------
try
    spm_vec(M.P) - spm_vec(M.pE);
catch
    M.P = M.pE;
end
try
    spm_vec(M.G) - spm_vec(M.gE);
catch
    M.G = M.gE;
end

% time-step
%--------------------------------------------------------------------------
try
    Y.dt;
catch
    Y.dt = 1;
end


% precision components Q
%--------------------------------------------------------------------------
try
    Q = Y.Q;
catch
    Q = spm_Ce(ns*ones(1,nr));
end
nh    = length(Q);          % number of precision components
ne    = length(Q{1});       % number of error terms
h     = zeros(nh,1);        % initialise hyperparameters
 

% confounds (if specified)
%--------------------------------------------------------------------------
try
    dgdu = kron(speye(nr,nr),Y.X0);
catch
    dgdu = sparse(ns*nr,0);
end

% hyperpriors - expectation
%--------------------------------------------------------------------------
try
    hE = M.hE;
catch
    hE = sparse(nh,1);
end

% hyperpriors - covariance
%--------------------------------------------------------------------------
try
    hP = inv(M.hC);
catch
    hP = eye(nh,nh);
end

% dimension reduction of parameter space
%--------------------------------------------------------------------------
Vp    = spm_svd(M.pC,exp(-16));
Vg    = spm_svd(M.gC,exp(-16));
np    = size(Vp,2);                   % number of parameters (f)
ng    = size(Vg,2);                   % number of parameters (g)
nu    = size(dgdu,2);                 % number of parameters (u)

% prior moments
%--------------------------------------------------------------------------
pE    = M.pE;
gE    = M.gE;
uE    = sparse(nu,1);
 
% second-order moments (in reduced space)
%--------------------------------------------------------------------------
pC    = Vp'*M.pC*Vp;
gC    = Vg'*M.gC*Vg;
uC    = speye(nu,nu)*exp(16);
ipC   = inv(pC);
igC   = inv(gC);
iuC   = inv(uC);
icC   = spm_cat(diag({igC,iuC}));
ibC   = spm_cat(diag({ipC icC}));

 
% initialise conditional density
%--------------------------------------------------------------------------
Ep    = M.P;
Eg    = M.G;
Eu    = uE;
 
% EM
%==========================================================================
C.F   = -Inf;
C.G   = -Inf;
tg    = 128;
tp    = 128;
dFdh  = zeros(nh,1);
dFdhh = zeros(nh,nh);
warning off

% Optimise p: parameters of f(x,u,p)
%==========================================================================
for k = 1:128
 
    
    % predicted hidden states (x) and dxdp
    %----------------------------------------------------------------------
    [dxdp x] = spm_diff(IS,Ep,M,U,1,{Vp});
    
       
    % Optimise g: parameters of G(g) (and confounds)
    %======================================================================
    for l = 1:8
        
        %  prediction yp = G(g)*x and errors
        %------------------------------------------------------------------
        [dGdg G] = spm_diff(GM,Eg,M,1,{Vg});
        yp       = x*G';
        
        ey       = spm_vec(y) - spm_vec(yp) - dgdu*Eu;
        ep       = Vp'*(spm_vec(Ep) - spm_vec(pE));
        eg       = Vg'*(spm_vec(Eg) - spm_vec(gE));
        eu       =      spm_vec(Eu) - spm_vec(uE);
        ec       = [eg; eu];

        % gradients
        %------------------------------------------------------------------
        for i = 1:np
            dgdp(:,i) = spm_vec(dxdp{i}*G');
        end
        for i = 1:ng
            dgdg(:,i) = spm_vec(x*dGdg{i}');
        end
        dgdc  = [dgdg dgdu];
        dgdb  = [dgdp dgdc];  

        % Optimise h: parameters of iS(h)
        %==================================================================
        for m = 1:8

            % precision
            %--------------------------------------------------------------
            iS    = speye(ne,ne)*1e-8;
            for i = 1:nh
                iS = iS + Q{i}*exp(h(i));
            end
            S     = inv(iS);
            dLdbb = dgdb'*iS*dgdb + ibC;
            Cb    = inv(dLdbb);

            % precision operators for M-Step
            %--------------------------------------------------------------
            for i = 1:nh
                P{i}  = Q{i}*exp(h(i));
                PS{i} = P{i}*S;
            end

            % derivatives: dLdh = dL/dh,...
            %--------------------------------------------------------------
            for i = 1:nh
                dFdh(i,1)      =  trace(PS{i})/2 - ey'*P{i}*ey/2 ...
                    -sum(sum(Cb.*(dgdb'*P{i}*dgdb)))/2;

                for j = i:nh
                    dFdhh(i,j) = -sum(sum(PS{i}.*PS{j}))/2;
                    dFdhh(j,i) =  dFdhh(i,j);
                end
            end

            % add hyperpriors
            %--------------------------------------------------------------
            dFdh  = dFdh  - hP*(h - hE);
            dFdhh = dFdhh - hP;

            % M-Step: update ReML estimate of h
            %--------------------------------------------------------------
            Ch    = inv(-dFdhh);
            dh    = Ch*dFdh;
            h     = h  + dh;

            % prevent overflow
            %--------------------------------------------------------------
            h     = max(h,-16);

            % convergence
            %--------------------------------------------------------------
            dF    = dFdh'*dh;
            if dF < 1e-2, break, end

        end
        
            % objective function: F(p) (= log-evidence - divergence)
        %==================================================================
        F = ...
        - ey'*iS*ey/2 ...
        - ep'*ipC*ep/2 ...
        - eg'*igC*eg/2 ...
        - eu'*iuC*eu/2 ...
        - ne*log(8*atan(1))/2 ...
        + spm_logdet(iS )/2 ...
        + spm_logdet(ipC)/2 ...
        + spm_logdet(igC)/2 ...
        + spm_logdet(iuC)/2 ...
        + spm_logdet(Cb )/2 ...
        + spm_logdet(Ch )/2;
    

        % if F has increased, update gradients and curvatures for E-Step
        %------------------------------------------------------------------
        if F > C.G

            % update gradients and curvature
            %--------------------------------------------------------------
            dFdc  =  dgdc'*iS*ey   - icC*ec;
            dFdcc = -dgdc'*iS*dgdc - icC;

            % accept current estimates
            %--------------------------------------------------------------
            C.Eg  = Eg;
            C.Eu  = Eu;
            C.h   = h;
            C.G   = F;

            % and decrease regularization
            %--------------------------------------------------------------
            tg    = tg*2;
            
        else

            % reset expansion point
            %--------------------------------------------------------------
            Eg    = C.Eg;
            Eu    = C.Eu;
            h     = C.h;

            % and increase regularization
            %--------------------------------------------------------------
            tg    = min(tg/2,1);

        end
        
        % E-Step: Conditional updates of g and u
        %------------------------------------------------------------------
        dc    = spm_dx(dFdcc,dFdc,{tg});
        dg    = dc(1:ng);
        du    = dc([1:nu]+ ng);
       
        Eg    = spm_unvec(spm_vec(Eg) + Vg*dg,Eg);
        Eu    = spm_unvec(spm_vec(Eu) + du,Eu);
       
        % convergence
        %------------------------------------------------------------------
        dF    = dFdc'*dc;
        if dF < 1e-2, break, end
        
    end
    
    % objective function: F(p) (= log-evidence - divergence)
    %======================================================================
    F = - ey'*iS*ey/2 ...
        - ep'*ipC*ep/2 ...
        - eg'*igC*eg/2 ...
        - eu'*iuC*eu/2 ...
        - ne*log(8*atan(1))/2 ...
        + spm_logdet(iS )/2 ...
        + spm_logdet(ipC)/2 ...
        + spm_logdet(igC)/2 ...
        + spm_logdet(iuC)/2 ...
        + spm_logdet(Cb )/2 ...
        + spm_logdet(Ch )/2;
 
 
    % if F has increased, update gradients and curvatures for E-Step
    %----------------------------------------------------------------------
    if F > C.F            
        
        % update gradients and curvature
        %------------------------------------------------------------------
        dFdp  =  dgdp'*iS*ey   - ipC*ep;
        dFdpp = -dgdp'*iS*dgdp - ipC;
        
        % accept current estimates
        %------------------------------------------------------------------
        C.Ep  = Ep;
        C.Eg  = Eg;
        C.Eu  = Eu;
        C.h   = h;
        C.F   = F;
 
        % and decrease regularization
        %------------------------------------------------------------------
        tp    = tp*2;
        str   = 'EM-Step(-)';
 
    else
 
        % reset expansion point
        %------------------------------------------------------------------
        Ep    = C.Ep;
        Eg    = C.Eg;
        Eu    = C.Eu;
        h     = C.h;
 
        % and increase regularization
        %------------------------------------------------------------------
        tp    = min(tp/2,1);
        str   = 'EM-Step(+)';
 
    end
 
    % Optimise p: parameters of f(x,u,p)
    %======================================================================
    dp    = spm_dx(dFdpp,dFdp,{tp});
    Ep    = spm_unvec(spm_vec(Ep) + Vp*dp,Ep);

    % graphics
    %----------------------------------------------------------------------
    try
        figure(Fsi)

        % subplot prediction
        %------------------------------------------------------------------
        subplot(2,1,1)
        plot([1:ns]*Y.dt,yp),                        hold on
        plot([1:ns]*Y.dt,yp + spm_unvec(ey,yp),':'), hold off
        xlabel('time')
        title(sprintf('%s: %i','E-Step',k))
        grid on

        % subplot parameters - f(P)
        %------------------------------------------------------------------
        subplot(2,2,3)
        bar(full(Vp*ep))
        xlabel('parameter f(x)')
        title('conditional [minus prior] expectation')
        grid on
        
        % subplot parameters - g(G)
        %------------------------------------------------------------------
        subplot(2,2,4)
        bar(full(Vg*eg))
        xlabel('parameter (g(x))')
        title('conditional [minus prior] expectation')
        grid on
        drawnow
    end

    % convergence
    %----------------------------------------------------------------------
    dF  = dFdp'*dp;
    fprintf('%-6s: %i %6s %e %6s %e\n',str,k,'F:',C.F,'dF:',full(dF))
    if k > 2 && dF < 1e-2, break, end

end
 
% outputs
%--------------------------------------------------------------------------
Cp     = Vp*Cb([1:np],[1:np])*Vp';
Cg     = Vg*Cb([1:ng] + np,[1:ng] + np)*Vg';
F      = C.F;
warning on




