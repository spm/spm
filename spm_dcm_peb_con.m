function [P,c,v] = spm_dcm_peb_con(PEB, C, threshold, doplot)
% Bayesian contrast of the parameters in a PEB model or BMA
% FORMAT [P,c,v] = spm_dcm_peb_con(PEB, C, threshold, doplot)
%
% Inputs:
%
% PEB       - estimated PEB model or BMA of PEB models
% C         - contrast vector or matrix (see below)
% threshold - (optional) test statistic [default: 0]
% doplot    - (optional) whether to plot results, default false
%
% Outputs:
% 
% P         - probability that the contrast value is larger than the
%             threshold
% c,v       - mean and variance of the contrast, e.g. the probability of a
%             difference between two connections
%
% Specifying the contrast vector / matrix:
%
% The contrast vector or matrix should be the same size as the parameters
% in the PEB model (if shorter, it will be padded with zeros). The 
% parameters in matrix PEB.Ep are of dimension: [connections x covariates].
% This matrix is vectorized in a BMA, by stacking the covariates on top of  
% one another. The elements of the contrast matrix define a linear mixture
% of the parameters. For example, for a PEB model with two connections and 
% three covariates, the following contrast compares the first and second 
% connections of covariate one:
%
% C = [1 0 0;
%     -1 0 0];
%
% And the following contrast compares the effects of the second and third 
% covariates on the first connection:
%
% C = [0 1 -1;
%      0 0  0];
% 
% Having defined the contrast matrix, call this function using:
%
% [P,c,v] = spm_dcm_peb_con(PEB, C, 0, true);
%__________________________________________________________________________

% Peter Zeidman
% Copyright (C) 2020-2022 Wellcome Centre for Human Neuroimaging

% Validate
if length(PEB) ~= 1
    error('Please provide a single PEB model, e.g. PEB(1)');
end
if nargin < 3, threshold = 0;  end
if nargin < 4, doplot = false; end

% Unpack PEB
Ep = spm_vec(PEB.Ep);
Cp = PEB.Cp;

% Variance vector -> covariance matrix
if length(Ep) > 1 && size(Cp,2) == 1
    Cp = diag(Cp);
end

% Vectorize contrast and zero pad
C   = C(:);
pad = zeros(length(Ep) - length(C),1);
C   = vertcat(C, pad);
if length(C) > length(Ep)
    error('The contrast vector is longer than parameter vector');
end

% Compute contrast
c    = C'*Ep;
v    = C'*Cp*C;

% Compute probability of a difference
P = 1 - spm_Ncdf(threshold,c,v);

% Unvectorize contrast for convenience
C = spm_unvec(C, PEB.Ep);

% Stop if not plotting
if ~doplot, return; end

% Create figure
Fgraph = spm_figure('GetWin','Graphics');
spm_figure('Clear',Fgraph);

% Compute normal density (code originally from spm_dcm_review)
x    = c + [-512:512]*sqrt(v)*6/512;
p    = full(1/sqrt(2*pi*v)*exp(-[x - c].^2/(2*v)));
sw   = warning('off','SPM:negativeVariance');
warning(sw);

% Plot density
figure(Fgraph)
subplot(2,1,1)
plot(x,p,[1 1]*threshold,[0 max(p)],'-.');
str  = sprintf('%s P(contrast > %0.2f) = %.1f%s','Posterior density',threshold,P*100,'%');
xlabel('contrast of parameter estimates')
ylabel('probability density')
set(gca,'FontSize',12);
title(str,'FontSize',16);

% Shade area under curve that's above threshold
i    = find(x >= threshold);
hold on
fill([x(i) fliplr(x(i))],[i*0 fliplr(p(i))],[1 1 1]*.8)
axis square, grid on
hold off

% Plot contrast matrix
subplot(2,1,2)
imagesc(C);
title('Contrast','FontSize',12);
if size(PEB.Ep > 1)
    try, set(gca,'XTick',1:length(PEB.Xnames),'XTickLabel',PEB.Xnames); end
    try, set(gca,'YTick',1:length(PEB.Pnames),'YTickLabel',PEB.Pnames); end
    xlabel('Covariate'); ylabel('Connection');
else
    ylabel('Parameter');
end
axis image
set(gca,'FontSize',12);
