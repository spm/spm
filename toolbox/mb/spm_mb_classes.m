function varargout = spm_mb_classes(varargin)
% Get tissue classes
% FORMAT [P,dat] = spm_mb_classes(dat,mu,sett)
% dat  - Data structure for a subject
% mu   - Warped template data
% sett - Settings
% P    - Updated tissue classes
%
% FORMAT [dat,P] = spm_mb_classes('update_cat',dat,mu,sett)
% FORMAT lab     = spm_mb_classes('get_labels',dat,K1)
% FORMAT cm      = spm_mb_classes('get_label_conf_matrix',cm_map,w,K1)
% FORMAT l       = spm_mb_classes('LSE',mu,ax)
% FORMAT mu      = spm_mb_classes('template_k1',mu,ax)
% FORMAT l       = spm_mb_classes('LSE1',mu,ax)
%__________________________________________________________________________
% Copyright (C) 2019-2020 Wellcome Centre for Human Neuroimaging

% $Id: spm_mb_classes.m 7936 2020-08-20 10:12:39Z mikael $

if isa(varargin{1},'char')
    [varargout{1:nargout}] = spm_subfun(localfunctions,varargin{:});
else
    [varargout{1:nargout}] = get_classes(varargin{:});
end
%==========================================================================

%==========================================================================
function [P,dat] = get_classes(dat,mu,sett)

if isfield(dat.model,'cat')
    % Categorical model
    [dat,P] = update_cat(dat,mu);
elseif isfield(dat.model,'gmm')
    % GMM model

    % Expand mu to include the background class and combine with labels if
    % required.
    mu  = template_k1(mu);
    lab = get_labels(dat,size(mu,4));
    if numel(lab)>1
        % Add labels to template
        mu = mu + lab;
    end
    clear lab

    [dat,P] = spm_mb_appearance('update',dat,mu,sett);
else
    error('This should not happen');
end
%==========================================================================

%==========================================================================
function [dat,P] = update_cat(dat,mu)
% Categorical model
P  = spm_mb_io('get_data',dat.model.cat.f);
sk = dat.samp;
P  = P(1:sk:end,1:sk:end,1:sk:end,:);

% Compute subject-specific categorical cross-entropy loss between
% segmentation and template
msk      = all(isfinite(P),4) & all(isfinite(mu),4);
tmp      =  sum(P.*mu,4) - LSE(mu,4);
dat.E(1) = -sum(tmp(msk(:)));
dat.nvox =  sum(msk(:));
%==========================================================================

%==========================================================================
function lab = get_labels(dat, K1)
if isempty(dat.lab), lab = 0; return; end

% Load labels
lab    = spm_mb_io('get_data', dat.lab.f);
sk     = dat.samp;
lab    = lab(1:sk(1):end, 1:sk(2):end, 1:sk(3):end);
dm     = [size(lab) 1 1];
lab    = round(lab(:));
cm_map = dat.lab.cm_map; % cell array that defines the confusion matrix
lab(~isfinite(lab) | lab<1 | lab>numel(cm_map)) = numel(cm_map) + 1; % Prevent crash

% Get confusion matrix that maps from label value to (log) probability value
cm = get_label_conf_matrix(cm_map, dat.lab.w, K1);

% Build "one-hot" representation using confusion matrix
lab = reshape(cm(lab,:), [dm(1:3) K1]);
%==========================================================================

%==========================================================================
function cm = get_label_conf_matrix(cm_map, w, K1)
% FORMAT CM = get_label_cm(cm_map, w, K1)
% cm_map - Defines the confusion matrix
% w      - Weighting probability
% K1     - Number of classes
% cm     - confusion matrix
%
% Build Rater confusion matrix for one subject.
% This matrix maps template classes to manually segmented classes.
% Manual labels often do not follow the same convention as the Template,
% and not all regions may be labelled. Therefore, a manual label may
% correspond to several Template classes and, conversely, one Template
% class may correspond to several manual labels.

% Parse function settings
w  = min(max(w, 1e-7), 1-1e-7);
L  = numel(cm_map);             % Number of labels
cm = zeros([L + 1 K1], 'single'); % Allocate confusion matrix (including unknown)
for l=1:L % Loop over labels
    ix            = false(1, K1);
    ix(cm_map{l}) = true;
    cm(l, ix)     = log(w/sum( ix));
    cm(l,~ix)     = log((1-w)/sum(~ix));
end
%==========================================================================

%==========================================================================
function l = LSE(mu,ax)
% log-sum-exp function
if nargin<2, ax = 4; end
mx = max(max(mu,[],ax),0);
l  = log(exp(-mx) + sum(exp(bsxfun(@minus,mu,mx)),ax)) + mx;
%==========================================================================

%==========================================================================
function mu = template_k1(mu,ax)
% Expand a template to include the implicit background class
if nargin<2, ax = 4; end
lse  = LSE(mu,ax);
mu   = cat(ax,bsxfun(@minus,mu,lse), -lse);
%==========================================================================

%==========================================================================
function l = LSE1(mu,ax)
% log-sum-exp (including final class) function
if nargin<2, ax = 4; end
mx = max(mu,[],ax);
l  = log(sum(exp(bsxfun(@minus,mu,mx)),ax)) + mx;
%==========================================================================

%==========================================================================

%==========================================================================

%==========================================================================

