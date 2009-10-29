function [Ep,Eg,Cp,Cg,S,F] = spm_nlsi_Nf(M,U,Y)
% Bayesian inversion of a nonlinear model of the form G(g)*F(u,p)
% FORMAT [Ep,Eg,Cp,Cg,S,F]= spm_nlsi_N(M,U,Y)
%
% Dynamical MIMO models
%__________________________________________________________________________
%
% M.IS - A prediction generating function name; this is usually an
%        integration scheme for states-sapce models of the form
%
%     M.f  - f(x,u,p,M) - state equation:  dxdt = f(x,u)
%     M.G  - G(x,M)     - linear observer: y    = G*x
%
% M.FS - function name f(y,M)   - feature selection
%        This [optional] function performs feature selection assuming the
%        generalized model y = FS(y,M) = FS(G*x,M) + X0*P0 + e
%
% M.x  - The expansion point for the states (i.e., the fixed point)
%
% M.P  - starting estimtes for model parameters [optional]
%
% M.pE - prior expectation  - of model parameters - f(x,u,p,M)
% M.pC - prior covariance   - of model parameters - f(x,u,p,M)
%
% M.gE - prior expectation  - of model parameters - G(x,M)
% M.gC - prior covariance   - of model parameters - G(x,M)
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
% If the free-energy starts to increase, a Levenberg-Marquardt scheme is
% invoked.  The M-Step estimates the precision components of e, in terms
% of [Re]ML point estimators of the log-precisions.
% An optional feature selection can be specified with parameters M.FS
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_nlsi_Nf.m 3517 2009-10-29 15:11:56Z guillaume $


if isfield(M,'fastEM') && M.fastEM
    % Check parameters structure...
    P0.T = [];
    P0.H = [];
    P0.S = [];
    P0.A = [];
    P0.B = [];
    P0.C = [];
    P0.D = [];
    P0.R = [];
    M.pE = orderfields(M.pE,P0);
    if ~isempty(M.P)
        M.P = orderfields(M.P,P0);
    end
else
    M.fastEM = 0;
end


% figure (unless disabled)
%--------------------------------------------------------------------------
try
    M.nograph;
catch
    Fsi = spm_figure('GetWin','SI');
end

% check integrator
%--------------------------------------------------------------------------
try
    IS = M.IS;
catch
    IS = 'spm_int_U';
end

% check observer has not been accidentally specified
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
        FS = inline([M.FS '(y,M)'],'y','M');

        % FS(y)
        %----------------------------------------------------------------------
    catch
        y  = feval(M.FS,Y.y);
        FS = inline([M.FS '(y)'],'y','M');

    end
else

    % y
    %----------------------------------------------------------------------
    y  = Y.y;
    FS = inline('y','y','M');
end

% size of data (usually samples x channels)
%--------------------------------------------------------------------------
if iscell(y)

    ns = size(y{1},1);

    % concatenate samples over cell, sensuring the same for prections
    %----------------------------------------------------------------------
    y  = spm_cat(y(:));
    IS = inline(['spm_cat(' IS '(P,M,U))'],'P','M','U');

else
    ns = size(y,1);
end
nr   = length(spm_vec(y))/ns;            % number of samples and responses
M.ns = ns;                               % store in M.ns for integrator

% initial states
%--------------------------------------------------------------------------
try
    M.x;
catch
    try
        M.n;
    catch
        M.n = 0;
    end
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
    iP    = spm_vec(M.P);
    pE    = spm_vec(M.pE);
    i     = find(diag(M.pC));
    pE(i) = iP(i);
    M.P   = spm_unvec(pE,M.pE);
    fprintf('\nParameter initialisation successful\n')
catch
    M.P = M.pE;
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
nt    = length(Q{1});       % number of time bins
nq    = nr*ns/nt;           % for compact Kronecker form of M-step
h     = zeros(nh,1);        % initialize hyperparameters


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
    hE  = M.hE;
catch
    hE  = sparse(nh,1) - 32;
end

% hyperpriors - covariance
%--------------------------------------------------------------------------
try
    ihC = inv(M.hC);
catch
    ihC = speye(nh,nh)/256;
end

% dimension reduction of parameter space
%--------------------------------------------------------------------------
Vp    = spm_svd(M.pC,exp(-32));
Vg    = spm_svd(M.gC,exp(-32));
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
sw    = warning('off','all');
pC    = Vp'*M.pC*Vp;
gC    = Vg'*M.gC*Vg;
uC    = speye(nu,nu)*exp(32);
ipC   = inv(pC);                      % p - state parameters
igC   = inv(gC);                      % g - observer parameters
iuC   = inv(uC);                      % u - fixed parameters
ibC   = spm_cat(spm_diag({ipC,igC,iuC})); % b - all parameters


% initialize conditional density
%--------------------------------------------------------------------------
Ep    = M.P;
Eg    = M.gE;
Eu    = spm_pinv(dgdu)*spm_vec(y);
warning(sw);

% EM
%==========================================================================
C.F   = -Inf;                                 % free-energy: f(x,u,p)
v     = 0;                                    % ascent rate: f(x,u,p)
dFdh  = zeros(nh,1);
dFdhh = zeros(nh,nh);
sw    = warning('off','all');

tic

% expansion point
%--------------------------------------------------------------------------
x0     = ones(size(y,1),1)*spm_vec(M.x)';


if ~M.fastEM
    disp('--- Nonlinear system identification: EM probabilistic inversion ----')
else
    disp('--- Nonlinear system identification: fast EM probabilistic inversion ----')
end


% Optimize p: parameters of f(x,u,p)
%==========================================================================
for ip = 1:64

    now = clock;

    % predicted hidden states (x) and dxdp
    %----------------------------------------------------------------------
    if ~M.fastEM
        [dxdp x] = spm_diff(IS,Ep,M,U,1,{Vp});
    else % fast EM
        M.Vp = Vp;
        [x,dxdp] = spm_gen_erp(Ep,M,U);
        x = spm_cat(x);
        nc = length(dxdp);
        dxdp2 = cell(nc,np);
        for c = 1:nc
            dxdp{c,1} = dxdp{c,1}*Vp;
            for i = 1:np
                dxdp2{c,i} = reshape(dxdp{c,1}(:,i),M.n,M.ns)';
            end
        end
        dxdp = spm_cat(dxdp2,1);
        clear dxdp2;
    end
    % check for dissipative dynamics
    %----------------------------------------------------------------------
    if all(isfinite(spm_vec(x)))
        gi   = 8;
    else
        gi   = 0;
    end


    % Optimize g: parameters of G(g)
    %======================================================================
    for ig = 1:gi

        % prediction yp = G(g)*x
        %------------------------------------------------------------------
        [dGdg G] = spm_diff(M.G,Eg,M,1,{Vg});
        yp       = FS((x - x0)*G',M);


        % Optimize g: parameters confounds
        %==================================================================

        % prediction errors
        %------------------------------------------------------------------
        ey    = spm_vec(y)  - spm_vec(yp) - dgdu*Eu;
        eu    = spm_vec(Eu) - spm_vec(uE);

        % update gradients and curvature
        %------------------------------------------------------------------
        dFdu  =  dgdu'*ey   - iuC*eu;
        dFduu = -dgdu'*dgdu - iuC;

        % Conditional updates of  u
        %------------------------------------------------------------------
        du    = spm_dx(dFduu,dFdu);
        Eu    = Eu + du;

        % recompute prediction errors
        %------------------------------------------------------------------
        ey    = spm_vec(y) - spm_vec(yp) - dgdu*Eu;
        ep    = Vp'*(spm_vec(Ep) - spm_vec(pE));
        eg    = Vg'*(spm_vec(Eg) - spm_vec(gE));
        eu    =      spm_vec(Eu) - spm_vec(uE);
        eb    = [ep; eg; eu];


        % gradients
        %------------------------------------------------------------------

        for i = 1:np
            dgdp(:,i) = spm_vec(FS(dxdp{i}*G',M));
        end

        try
            for i = 1:ng
                dgdg(:,i) = spm_vec(FS((x - x0)*dGdg{i}',M));
            end
        catch
            dgdg = FS((x - x0)*dGdg,M);
        end
        dgdb  = [dgdp dgdg dgdu];

        % Optimize F(h): parameters of iS(h)
        %==================================================================
        for ih = 1:8

            % precision
            %--------------------------------------------------------------
            iS    = speye(nt,nt)*exp(-32);
            for i = 1:nh
                iS = iS + Q{i}*exp(h(i));
            end
            S     = inv(iS);
            iS    = kron(speye(nq),iS);
            dFdbb = dgdb'*iS*dgdb + ibC;
            Cb    = inv(dFdbb);

            % precision operators for M-Step
            %--------------------------------------------------------------
            for i = 1:nh
                P{i}  = Q{i}*exp(h(i));
                PS{i} = P{i}*S;
                P{i}  = kron(speye(nq),P{i});
            end

            % derivatives: dLdh = dL/dh,...
            %--------------------------------------------------------------
            for i = 1:nh
                dFdh(i,1)      =  trace(PS{i})*nq/2 - ey'*P{i}*ey/2 ...
                    -sum(sum(Cb.*(dgdb'*P{i}*dgdb)))/2;
                for j = i:nh
                    dFdhh(i,j) = -sum(sum(PS{i}.*PS{j}))*nq/2;
                    dFdhh(j,i) =  dFdhh(i,j);
                end
            end

            % add hyperpriors
            %--------------------------------------------------------------
            eh    = h     - hE;
            dFdh  = dFdh  - ihC*eh;
            dFdhh = dFdhh - ihC;
            Ch    = inv(-dFdhh);

            % M-Step: update ReML estimate of h
            %--------------------------------------------------------------
            dh    = spm_dx(dFdhh,dFdh);
            h     = h + dh;

            % prevent overflow
            %--------------------------------------------------------------
            h     = min(max(h,-16),16);

            % convergence
            %--------------------------------------------------------------
            if dFdh'*dh < 1e-1, break, end

        end

        % optimise F(g,u)
        %==================================================================

        % update gradients and curvature
        %------------------------------------------------------------------
        dFdg  =  dgdg'*iS*ey   - igC*eg;
        dFdgg = -dgdg'*iS*dgdg - igC;

        % E-step: Conditional updates of g and u
        %------------------------------------------------------------------
        dg    = spm_dx(dFdgg,dFdg);
        Eg    = spm_unvec(spm_vec(Eg) + Vg*dg,Eg);

        % convergence
        %------------------------------------------------------------------
        dG    = dFdg'*dg;
        if ig > 1 && dG < 1e-1, break, end

    end

    % optimise objective function: F(p) = log-evidence - divergence
    %======================================================================
    F = ...
        - ey'*iS*ey/2 ...
        - ep'*ipC*ep/2 ...
        - eg'*igC*eg/2 ...
        - eu'*iuC*eu/2 ...
        - eh'*ihC*eh/2 ...
        - ns*nr*log(8*atan(1))/2 ...
        - nq*spm_logdet(S)/2 ...
        + spm_logdet(ibC*Cb)/2 ...
        + spm_logdet(ihC*Ch)/2;

    % record increases and reference log-evidence for reporting
    %----------------------------------------------------------------------
    try
        F0; fprintf(' actual: %.3e',full(F - C.F))
        elapsedTime = etime(clock,now);
        fprintf(['  (took ',num2str(elapsedTime),' sec)\n'])
    catch
        F0 = F;
    end

    % if F has increased, update gradients and curvatures for E-Step
    %----------------------------------------------------------------------
    if F > C.F

        % update gradients and curvature
        %------------------------------------------------------------------
        dFdp  =  dgdp'*iS*ey   - ipC*ep;
        dFdpp = -dgdp'*iS*dgdp - ipC;

        % decrease regularization
        %------------------------------------------------------------------
        v     = min(v + 1,8);
        str   = 'EM(+)';

        % accept current estimates
        %------------------------------------------------------------------
        C.Cb  = Cb;                               % conditional covaraince
        C.Ep  = Ep;                               % and expectations
        C.Eg  = Eg;
        C.Eu  = Eu;
        C.h   = h;
        C.F   = F;

    else

        % reset expansion point
        %------------------------------------------------------------------
        Cb    = C.Cb;                             % conditional covaraince
        Ep    = C.Ep;                             % and expectations
        Eg    = C.Eg;
        Eu    = C.Eu;
        h     = C.h;

        % and increase regularization
        %------------------------------------------------------------------
        v     = min(v - 1,0);
        str   = 'EM(-)';

    end

    % Optimize p: parameters of f(x,u,p)
    %======================================================================
    dp    = spm_dx(dFdpp,dFdp,{v});
    Ep    = spm_unvec(spm_vec(Ep) + Vp*dp,Ep);



    % subplot times
    %----------------------------------------------------------------------
    if length(Y.pst) == size(yp,1)
        yt = Y.pst;
    else
        yt = [1:size(yp,1)]*Y.dt*1000;
    end

    % graphics
    %----------------------------------------------------------------------
    try

        % subplot prediction
        %------------------------------------------------------------------
        figure(Fsi)

        subplot(3,1,1)
        plot(yt,x)
        xlabel('time (ms)')
        set(gca,'XLim',[yt(1) yt(end)])
        title(sprintf('%s: %i','E-Step: hidden states',ip))
        grid on

        subplot(3,1,2)
        plot(yt,yp),                        hold on
        plot(yt,yp + spm_unvec(ey,yp),':'), hold off
        xlabel('time (ms)')
        set(gca,'XLim',[yt(1) yt(end)])
        title('E-Step: response amd prediction')
        grid on

        % subplot parameters - f(P)
        %------------------------------------------------------------------
        subplot(3,2,5)
        bar(full(Vp*ep))
        xlabel('parameter f(x)')
        title('conditional [minus prior] expectation')
        grid on

        % subplot parameters - g(G)
        %------------------------------------------------------------------
        subplot(3,2,6)
        bar(full(Vg*eg))
        xlabel('parameter (g(x))')
        title('conditional [minus prior] expectation')
        grid on
        drawnow

    end

    % convergence
    %----------------------------------------------------------------------
    dF  = dFdp'*dp;
    ig  = max([0 ig]);
    fprintf('%-6s: %-2i (%i,%i) %4s %-6.3e %6s %6.3e ',str,ip,ig,ih,'F:',full(C.F - F0),'dF predicted:',full(dF))
    if ip > 2 && dF < 1e-2
        fprintf(' convergence\n')
        break
    end
end

% outputs
%--------------------------------------------------------------------------
Cp     = Vp*Cb([1:np],     [1:np]     )*Vp';
Cg     = Vg*Cb([1:ng] + np,[1:ng] + np)*Vg';
F      = C.F;
fprintf('\n')
toc
warning(sw);
return




% notes
%==========================================================================
Up    = spm_svd(dFdpp,0);
dFdp  = Up'*dFdp;
dFdpp = Up'*dFdpp*Up;
Vp    = Vp*Up;
