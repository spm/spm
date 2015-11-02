function [Ep,Cp,F] = spm_nlsi_Newton(M,U,Y)
% Variational Lapalce for nonlinear models - Newton's method
% FORMAT [Ep,Cp,F] = spm_nlsi_Newton(M,U,Y)
%
% Eplicit log-likihood model
%__________________________________________________________________________
%
% M.L - log likelihood function @(P,M,U,Y)
%       P  - free parameters
%       M  - model
%
% M.P  - starting estimates for model parameters [optional]
% M.pE - prior expectation      - E{P}   of model parameters
% M.pC - prior covariance       - Cov{P} of model parameters
%
% U  - inputs or causes
% Y  - output or response
%
% Parameter estimates
%--------------------------------------------------------------------------
% Ep  - (p x 1)         conditional expectation    E{P|y}
% Cp  - (p x p)         conditional covariance     Cov{P|y}
%
% log evidence
%--------------------------------------------------------------------------
% F   - [-ve] free energy F = log evidence = p(Y|pE,pC) = p(y|m)
%
%__________________________________________________________________________
% Returns the moments of the posterior p.d.f. of the parameters of a
% nonlinear model with a log likelihood function L(P,M,U,Y).
%
% Priors on the free parameters P are specified in terms of expectation pE
% and covariance pC. This Variational Laplace scheme uses an explicit
% (numerical) curvature to implement a gradient ascent on variational free
% energy using Newton's method. An example of its application is provided at
% the end of this routine using a simple general linear model. This example
% eschews the mean field approximation aassociated with standard
% inversions.
%
% For generic aspects of the scheme see:
%
% Friston K, Mattout J, Trujillo-Barreto N, Ashburner J, Penny W.
% Variational free energy and the Laplace approximation.
% NeuroImage. 2007 Jan 1;34(1):220-34.
%__________________________________________________________________________
% Copyright (C) 2001-2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_nlsi_Newton.m 6587 2015-11-02 10:29:49Z karl $

% options
%--------------------------------------------------------------------------
try, M.nograph; catch, M.nograph = 0;   end
try, M.noprint; catch, M.noprint = 0;   end
try, M.Nmax;    catch, M.Nmax    = 128; end

% converted to function handle
%--------------------------------------------------------------------------
L   = spm_funcheck(M.L);

% initial parameters
%--------------------------------------------------------------------------
try 
    M.P; fprintf('\nParameter initialisation successful\n')
catch
    M.P = M.pE;
end

% prior moments (assume uninformative priors if not specifed)
%--------------------------------------------------------------------------
pE     = M.pE;
try
    pC = M.pC;
catch
    np = spm_length(M.pE);
    pC = speye(np,np)*exp(16);
end

% unpack covariance
%--------------------------------------------------------------------------
if isstruct(pC);
    pC = spm_diag(spm_vec(pC));
end

% dimension reduction of parameter space
%--------------------------------------------------------------------------
V     = spm_svd(pC,0);

% second-order moments (in reduced space)
%--------------------------------------------------------------------------
pC    = V'*pC*V;
ipC   = inv(pC);

% initialize conditional density
%--------------------------------------------------------------------------
p     = V'*(spm_vec(M.P) - spm_vec(M.pE));
Ep    = spm_unvec(spm_vec(pE) + V*p,pE);

% figure (unless disabled)
%--------------------------------------------------------------------------
if ~M.nograph, Fsi = spm_figure('GetWin','SI'); clf, end


% Wariational Laplace
%==========================================================================
criterion = [0 0 0 0];
C.F   = -Inf;                                    % free energy
v     = -4;                                      % log ascent rate
for k = 1:M.Nmax
    
    % time
    %----------------------------------------------------------------------
    tStart = tic;
    
    % Log-likelihood  f, gradients; dfdp and curvature dfdpp
    %======================================================================
    [dfdpp,dfdp,f] = spm_diff(L,Ep,M,U,Y,[1 1],{V});
    dfdp           = dfdp';    
    dfdpp          = full(spm_cat(dfdpp'));
    
    % enure prior bounds on curvature
    %----------------------------------------------------------------------
    [E,D] = eig(dfdpp);
    D     = diag(D);
    dfdpp = E*diag(D.*(D < 0))*E';
    
    % condiitonal covariance
    %----------------------------------------------------------------------
    Cp    = inv(ipC - dfdpp);
    
    % Fre  energy: F(p) = log evidence - divergence
    %======================================================================
    F     = f - p'*ipC*p/2 + spm_logdet(ipC*Cp)/2;
    G(k)  = F;
 
    % record increases and reference log-evidence for reporting
    %----------------------------------------------------------------------
    if k > 1
        if ~M.noprint
            fprintf(' actual: %.3e (%.2f sec)\n',full(F - C.F),toc(tStart))
        end
    else
        F0 = F;
    end
    
    % if F has increased, update gradients and curvatures for E-Step
    %----------------------------------------------------------------------
    if F > C.F || k < 8
        
        % accept current estimates
        %------------------------------------------------------------------
        C.p   = p;
        C.F   = F;
        C.Cp  = Cp;
        
        % E-Step: Conditional update of gradients and curvature
        %------------------------------------------------------------------
        dFdp  = dfdp  - ipC*p;
        dFdpp = dfdpp - ipC;
        
        % decrease regularization
        %------------------------------------------------------------------
        v     = min(v + 1/2,4);
        str   = 'EM:(+)';
        
    else
        
        % reset expansion point
        %------------------------------------------------------------------
        p     = C.p;
        
        % and increase regularization
        %------------------------------------------------------------------
        v     = min(v - 2,-4);
        str   = 'EM:(-)';
        
    end
    
    % E-Step: update
    %======================================================================
    dp    = spm_dx(dFdpp,dFdp,{v});
    p     = p + dp;
    Ep    = spm_unvec(spm_vec(pE) + V*p,pE);
    
    % Graphics
    %======================================================================
    if exist('Fsi', 'var')
        
        spm_figure('Select', Fsi)
        
        % trajectory in parameter space
        %------------------------------------------------------------------
        subplot(2,2,1)
        plot(0,0,'r.','MarkerSize',32), hold on
        col = [exp(-k/4) exp(-k) 1];
        try
            plot(V(1,:)*p,V(2,:)*p,'.','MarkerSize',32,'Color',col), hold on
            xlabel('1st parameter')
            ylabel('2nd parameter')
        catch
            plot(k,V(1,:)*p,'.','MarkerSize',32,'Color',col), hold on
            xlabel('Iteration')
            ylabel('1st parameter')
        end
        title('Trajectory','FontSize',16)
        grid on, axis square
        
        % trajectory in parameter space
        %------------------------------------------------------------------
        subplot(2,2,2)
        bar(full(G - F0),'c')
        xlabel('Iteration')
        ylabel('Log-evidence')
        title('Free energy','FontSize',16)
        grid on, axis square
        
        % subplot parameters
        %--------------------------------------------------------------
        subplot(2,2,3)
        bar(full(spm_vec(pE) + V*p))
        xlabel('Parameter')
        tstr = 'Conditional expectation';
        title(tstr,'FontSize',16)
        grid on, axis square
        
        % subplot parameters (eigenmodes)
        %------------------------------------------------------------------
        subplot(2,2,4)
        spm_plot_ci(p,Cp)
        xlabel('Parameter (eigenmodes)')
        title('Posterior deviations','FontSize',16)
        grid on, axis square
        drawnow
        
    end
    
    % convergence
    %----------------------------------------------------------------------
    dF  = dFdp'*dp;
    if ~M.noprint
        fprintf('%-6s: %i %6s %-6.3e %6s %.3e ',str,k,'F:',full(C.F - F0),'dF predicted:',full(dF))
    end
    criterion = [(dF < 1e-1) criterion(1:end - 1)];
    if all(criterion)
        if ~M.noprint
            fprintf(' convergence\n')
        end
        break
    end
    
end

if exist('Fsi', 'var')
    spm_figure('Focus', Fsi)
end

% outputs
%--------------------------------------------------------------------------
Ep     = spm_unvec(spm_vec(pE) + V*C.p,pE);
Cp     = V*C.Cp*V';
F      = C.F;


return


% NB: notes - illustrative application (a simple linear model)
%==========================================================================

% parameters P and design matrix U
%--------------------------------------------------------------------------
U      = randn(32,2);                          % design matrix
P.beta = [4;2];                                % parameters of GLM
P.pi   = 2;                                    % log precision

% generate data
%--------------------------------------------------------------------------
Y      = U*P.beta + exp(-P.pi/2)*randn(32,1);

% model specification with log-likelihood function M.L
%--------------------------------------------------------------------------
M.L    = @(P,M,U,Y) sum(log( spm_Npdf(Y, U*P.beta, exp(-P.pi)) ));
M.pE   = spm_zeros(P);                           % prior means (parameters)
M.pC   = eye(spm_length(P));                     % prior variance (parameters)

% Variational Laplace
%--------------------------------------------------------------------------
[Ep,Cp,F] = spm_nlsi_Newton(M,U,Y);

% overlay true values on confidence intervals
%--------------------------------------------------------------------------
subplot(2,2,4),hold on
bar(spm_vec(P),1/4)



