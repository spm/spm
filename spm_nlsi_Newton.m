function [Ep,Cp,F] = spm_nlsi_Newton(M,U,Y)
% Bayesian inversion of nonlinear models - Newton's method
% FORMAT [Ep,Cp,F] = spm_nlsi_Newton(M,U,Y)
%
% Eplicit log-likihood models
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
% Y  - outputs or response
%
% Parameter estimates
%--------------------------------------------------------------------------
% Ep  - (p x 1)         conditional expectation    E{P|y}
% Cp  - (p x p)         conditional covariance     Cov{P|y}
%
% log evidence
%--------------------------------------------------------------------------
% F   - [-ve] free energy F = log evidence = p(y|f,g,pE,pC) = p(y|m)
%
%__________________________________________________________________________
% Returns the moments of the posterior p.d.f. of the parameters of a
% nonlinear model specified by L(P,M,U).
%
% Priors on the free parameters P are specified in terms of expectation pE
% and covariance pC.
%
% For generic aspects of the scheme see:
%
% Friston K, Mattout J, Trujillo-Barreto N, Ashburner J, Penny W.
% Variational free energy and the Laplace approximation.
% NeuroImage. 2007 Jan 1;34(1):220-34.
%
% This scheme handels complex data along the lines originally described in:
%
% Sehpard RJ, Lordan BP, and Grant EH.
% Least squares analysis of complex data with applications to permittivity
% measurements.
% J. Phys. D. Appl. Phys 1970 3:1759-1764.
%
%__________________________________________________________________________
% Copyright (C) 2001-2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_nlsi_Newton.m 6586 2015-10-31 12:02:44Z karl $

% options
%--------------------------------------------------------------------------
try, M.nograph; catch, M.nograph = 0;   end
try, M.noprint; catch, M.noprint = 0;   end
try, M.Nmax;    catch, M.Nmax    = 128; end

% converted to function handle
%--------------------------------------------------------------------------
L   = spm_funcheck(M.L);

% size of data (samples x response component x response component ...)
%--------------------------------------------------------------------------
ny   = spm_length(Y);

% input
%--------------------------------------------------------------------------
try, U; catch, U = []; end

% initial parameters
%--------------------------------------------------------------------------
try 
    M.P; fprintf('\nParameter initialisation successful\n')
catch
    M.P = M.pE;
end


% prior moments (assume uninformative priors if not specifed)
%--------------------------------------------------------------------------
pE       = M.pE;
try
    pC   = M.pC;
catch
    np   = spm_length(M.pE);
    pC   = speye(np,np)*exp(16);
end



% unpack covariance
%--------------------------------------------------------------------------
if isstruct(pC);
    pC = spm_diag(spm_vec(pC));
end

% dimension reduction of parameter space
%--------------------------------------------------------------------------
V     = spm_svd(pC,0);
np    = size(V,2);                    % number of parameters (effective)


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
if ~M.nograph
      
    % trajectory in parameter space
    %----------------------------------------------------------------------
    Fsi = spm_figure('GetWin','SI');clf
    subplot(2,2,1)
    plot(p(1),p(2),'.','MarkerSize',16), hold on
    xlabel('1st parameter (eigenmode)')
    ylabel('2nd parameter (eigenmode)')
    title('Trajectory','FontSize',16)
    grid on, axis square
    
end


% EM
%==========================================================================
criterion = [0 0 0 0];

C.F   = -Inf;                                   % free energy
v     = 0;                                      % log ascent rate
for k = 1:M.Nmax
    
    % time
    %----------------------------------------------------------------------
    tStart = tic;
    
    % Log-likelihood  f, gradients; dfdp and curvature dfdpp
    %======================================================================
    [dfdpp,dfdp,f] = spm_diff(L,Ep,M,U,Y,[1 1],{V});
    dfdp           = dfdp(:);    
    dfdpp          = full(spm_cat(dfdpp(:)));
    
    % enure prior bounds on curvature
    %----------------------------------------------------------------------
    [E,D]  = eig(dfdpp);
    dfdpp  = E*D.*(D < 0)*E;
    
    % condiitonal covariance
    %----------------------------------------------------------------------
    Cp = inv(ipC - E*D*E);
    
    % Fre  energy: F(p) = log evidence - divergence
    %======================================================================
    F     = f - p'*ipC*p/2 + spm_logdet(ipC*Cp)/2;
    G(k)  = F;
 
    % record increases and reference log-evidence for reporting
    %----------------------------------------------------------------------
    if exist('F0','var')
        if ~M.noprint
            fprintf(' actual: %.3e (%.2f sec)\n',full(F - C.F),toc(tStart))
        end
    else
        F0 = F;
    end
    
    % if F has increased, update gradients and curvatures for E-Step
    %----------------------------------------------------------------------
    if F > C.F || k < 4
        
        % accept current estimates
        %------------------------------------------------------------------
        C.p   = p;
        C.F   = F;
        C.Cp  = Cp;
        
        % E-Step: Conditional update of gradients and curvature
        %------------------------------------------------------------------
        dFdp  = dfdp - ipC*p;
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
        plot(p(1),p(2),'.','MarkerSize',16), hold on
        xlabel('1st parameter (eigenmode)')
        ylabel('2nd parameter (eigenmode)')
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
        
        % plot real or complex predictions
        %------------------------------------------------------------------
        subplot(2,2,3)
        spm_plot_ci(p,Cp)
        xlabel('Parameter (eigenmode)')
        title('Posterior deviations','FontSize',16)
        grid on, axis square
            
        % subplot parameters
        %--------------------------------------------------------------
        subplot(2,2,4)
        bar(full(spm_vec(pE) + V*p))
        xlabel('Parameter')
        tstr = 'Conditional expectation';
        title(tstr,'FontSize',16)
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

