function [p] = spm_mvb_p(MVB,k)
% Classical p-value for MVB using null distribution of log-odds ratio
% FORMAT [p] = spm_mvb_p(MVB,k)
%
% MVB - Multivariate Bayes structure
% k   - number of samples > 20
%
% p   - p-value: of (relative) F using an empirical null distribution
%
% spm_mvb_p evaluates an empirical null distribution for the (fee-energy)
% difference in log-evidences (the log odds ratio) by phase-shuffling the
% target vector and repeating the greedy search. It adds the p-value as a
% field (p_value) to MVB.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_mvb_p.m 3656 2009-12-23 20:17:30Z karl $
 
%-number of samples
%--------------------------------------------------------------------------
try
    k;
catch
    str   = 'number of samples for null distribution';
    k     = spm_input(str,'!+1','b',{'20','50','100'},[20 50 100]);
end
 
%-Get figure handles and set title
%--------------------------------------------------------------------------
Fmvb = spm_figure('GetWin','MVB');
spm_clf(Fmvb);
 
% get MVB results
%==========================================================================
try
    MVB;
catch
    mvb  = spm_select(1,'mat','please select models',[],pwd,'MVB_*');
    MVB  = load(mvb(1,:));
    MVB  = MVB.MVB;
end
 
 
% whiten target and predictor (X) variables (Y) (i.e., remove correlations)
%--------------------------------------------------------------------------
K     = MVB.K;
X     = K*MVB.X;
Y     = K*MVB.Y;
X0    = K*MVB.X0;
U     = MVB.M.U;
  
% create orthonormal projection to remove confounds
%--------------------------------------------------------------------------
Ns    = length(X);
X0    = spm_svd(X0);
R     = speye(Ns) - X0*X0';
R     = spm_svd(R);
Y     = R'*Y;
 
% F value (difference in log-evidence or log odds ratio)
%--------------------------------------------------------------------------
F     = MVB.M.F;
F     = max(F) - F(1);
 
% Randomisation testing
%==========================================================================
for i = 1:k
 
    % randomise target vector (using phase-shuffling if V ~= I)
    %----------------------------------------------------------------------
    if MVB.V(1,2)
        X0 = R'*spm_phase_shuffle(X);
    else
        X0 = R'*X(randperm(Ns),:);
    end
 
    % Optimise mapping
    %======================================================================
    M     = spm_mvb(X0,Y,[],U,[],MVB.Ni,MVB.sg);
 
    % record F and compute p-value
    %----------------------------------------------------------------------
    F0(i) = max(M.F) - M.F(1);
    p     = 1 - sum(F0 < F)/(1 + i);
    
    % plot
    %----------------------------------------------------------------------
    subplot(2,1,1)
    q     = ceil(i/4);
    hist(F0,q); hold on
    plot([F F],[0 q],'r'); hold off;
    str = sprintf('p < %0.4f (%i samples)',p,i);
    title({[MVB.name ' (' MVB.contrast ')']; str},'FontSize',16)
    xlabel('log odds ratio')
    ylabel('null distribution')
    axis square
    drawnow
    
end
 
 
% display and assign in base memory
%--------------------------------------------------------------------------
str = sprintf('randomisation p-value = %.4f',p);
xlabel({'log odds ratio';str},'FontSize',16)
disp(['Thank you; ' str])
MVB.p_value = p;
 
% save results
%--------------------------------------------------------------------------
try
   save(MVB.name,'MVB');
end
assignin('base','MVB',MVB)
