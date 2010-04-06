function spm_mvb_cvk_display(MVB)
% model display for MVB with cross-validation
% FORMAT spm_mvb_cvk_display(MVB)
% MVB  - multivariate Bayes structure, select one if not provided
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Christophe Phillips
% $Id: spm_mvb_cvk_display.m 3806 2010-04-06 14:42:32Z ged $

if nargin<1
    load(spm_select(1,'^MVB.*\.mat','Select MVB to display'))
end
if ~isfield(MVB,'cvk')
    error(['No crossvalidation data available. ' ...
        'Select another file or perform the crossvalidation']);
end

%-Get figure handles and set title
%--------------------------------------------------------------------------
Fmvb = spm_figure('GetWin','MVB');
spm_clf(Fmvb);
 
% get stuff in place for display
%--------------------------------------------------------------------------
K     = MVB.K;
X     = K*MVB.X;
X0    = orth(K*MVB.X0);
R     = speye(length(X)) - X0*X0';
R     = orth(R);
pX     = R*R'*X;
 
% plot validation
%--------------------------------------------------------------------------
subplot(2,2,1)
s      = 1:length(pX);
plot(s,pX,s,MVB.cvk.qX,'-.')
xlabel('sample')
ylabel('response (adjusted)')
title('cross-validation')
axis square
 
subplot(2,2,2)
plot(pX,MVB.cvk.qX,'.')
xlabel('true')
ylabel('predicted')
title(sprintf('p-value (parametric) = %.5f',MVB.p_value))
axis square
abc = axis;
hold on
plot([max(abc([1 3])) min(abc([2 4]))],[max(abc([1 3])) min(abc([2 4]))],'k')

% plot feature weights
%--------------------------------------------------------------------------
subplot(2,2,3)
imagesc(corrcoef(MVB.cvk.qE))
colorbar
caxis([0 1])
xlabel('bipartition (k)')
title({'correlations among';'k-fold feature weights'})
axis square
 
subplot(2,2,4)
if ~isempty(MVB.XYZ) && ~isempty(MVB.VOX)
    if isfield(MVB.cvk, 'P')
        spm_mip(prod(MVB.cvk.P,2),MVB.XYZ(1:3,:),MVB.VOX)
        title({[MVB.name ' (' MVB.contrast ')'];'prod( P(|weights| > 0) )'})
    else % reproduce plot from original spm_mvb_cvk2
        qe = mean(MVB.cvk.qE,2);
        qe = qe.*(qe > 0);
        spm_mip(qe,MVB.XYZ(1:3,:),MVB.VOX)
        title({[MVB.name ' (' MVB.contrast ')'];'mean (positive) weights'})
    end        
else
    % (Allows MVB for non-spatial data, excluding smooth or compact priors)
    title('No spatial info in MVB')
end
axis square
 
fprintf('\np-value = %.4f; classification: %.1f%%; R-squared %.1f%%\n', ...
    MVB.p_value,MVB.percent,MVB.R2)
