function spm_mvb_cvk_display(MVB)
% model display for MVB with cross-validation
% FORMAT spm_mvb_cvk_display(MVB)
% MVB  - multivariate Bayes structure, select one if not provided
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Christophe Phillips
% $Id: spm_mvb_cvk_display.m 2356 2008-10-17 15:08:19Z christophe $

if nargin<1
    load(spm_select(1,'^MVB.*\.mat','Select MVB to display'))
end
if ~isfield(MVB,'cvk')
    error('No crossvalidation data available. Select another file or perform the crossvalidation');
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
xlabel('biparititon (k)')
title({'correlations among';'k-fold feature weights'})
axis square
 
subplot(2,2,4)
spm_mip(prod(MVB.cvk.P,2),MVB.XYZ(1:3,:),MVB.VOX)
title({[MVB.name ' (' MVB.contrast ')'];'prod( P(|weights| > 0) )'})
axis square
 
fprintf('\np-value = %.4f; classification: %.1f%s; R-squared %.1f%s\n', ...
            MVB.p_value,MVB.percent,'%',MVB.R2,'%')