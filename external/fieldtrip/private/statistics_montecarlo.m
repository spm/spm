function [stat, cfg] = statistics_montecarlo(cfg, dat, design)

% STATISTICS_MONTECARLO performs a nonparametric statistical test by calculating 
% Monte-Carlo estimates of the significance probabilities and/or critical values from the 
% permutation distribution. This function should not be called
% directly, instead you should call the function that is associated with
% the type of data on which you want to perform the test.
%
% Use as:
%   stat = timelockstatistics(cfg, data1, data2, data3, ...)
%   stat = freqstatistics    (cfg, data1, data2, data3, ...)
%   stat = sourcestatistics  (cfg, data1, data2, data3, ...)
% where the data is obtained from TIMELOCKANALYSIS, FREQANALYSIS
% or SOURCEANALYSIS respectively, or from TIMELOCKGRANDAVERAGE,
% FREQGRANDAVERAGE or SOURCEGRANDAVERAGE respectively.
%
% The configuration can contain
%   cfg.statistic        = string, statistic to compute for each sample or voxel (see below)
%   cfg.design           = design matrix
%   cfg.numrandomization = number of randomizations, can be 'all'
%   cfg.correctm         = apply multiple-comparison correction, 'no', 'max', cluster', 'bonferoni', 'holms', 'fdr' (default = 'no')
%   cfg.alpha            = critical value for rejecting the null-hypothesis (default = 0.05)
%   cfg.tail             = -1, 1 or 0 (default = 0)
%   cfg.ivar             = number or list with indices, independent variable(s)
%   cfg.uvar             = number or list with indices, unit variable(s)
%   cfg.wvar             = number or list with indices, within-cell variable(s)
%   cfg.cvar             = number or list with indices, control variable(s)
%   cfg.feedback         = 'gui', 'text', 'textbar' or 'no' (default = 'text')
%
% If you use a cluster-based statistic, you can specify the following
% options that determine how the single-sample or single-voxel
% statistics will be thresholded and combined into one statistical
% value per cluster.
%   cfg.clusterstatistic = how to combine the single samples that belong to a cluster, 'maxsum', 'maxsize', 'wcm' (default = 'maxsum')
%   cfg.clusterthreshold = method for single-sample threshold, 'parametric', 'nonparametric_individual', 'nonparametric_common' (default = 'parametric')
%   cfg.clusteralpha     = for either parametric or nonparametric thresholding (default = 0.05)
%   cfg.clustercritval   = for parametric thresholding (default is determined by the statfun)
%   cfg.clustertail      = -1, 1 or 0 (default = 0)
%
% To include the channel dimension for clustering, you should specify
%   cfg.neighbours       = structure with the neighbours of each channel, see NEIGHBOURHOODSELECTION
% If the neighbourhood structure is empty, clustering will only be done
% in frequency and time (if available) and not over neighbouring channels.
%
% The statistic that is computed for each sample in each random reshuffling 
% of the data is specified as
%   cfg.statistic       = 'indepsamplesT'     independent samples T-statistic,
%                         'indepsamplesF'     independent samples F-statistic,
%                         'indepsamplesregrT' independent samples regression coefficient T-statistic,
%                         'indepsamplesZcoh'  independent samples Z-statistic for coherence,
%                         'depsamplesT'       dependent samples T-statistic,
%                         'depsamplesF'       dependent samples F-statistic,
%                         'depsamplesregrT'   dependent samples regression coefficient T-statistic,
%                         'actvsblT'          activation versus baseline T-statistic.
%
% You can also use a custom statistic of your choise that is sensitive
% to the expected effect in the data. You can implement the statistic
% in a "statfun" that will be called for each randomization. The
% requirements on a custom statistical function is that the function
% is called statfun_xxx, and that the function returns a structure
% with a "stat" field containing the single sample statistical values.
% Check the private functions statfun_xxx (e.g.  with xxx=tstat) for
% the correct format of the input and output.
%
% See also TIMELOCKSTATISTICS, FREQSTATISTICS, SOURCESTATISTICS

% Undocumented local options:
%   cfg.resampling       permutation, bootstrap
%   cfg.computecritval   yes|no, for the statfun
%   cfg.computestat      yes|no, for the statfun
%   cfg.computeprob      yes|no, for the statfun
%   cfg.voxelstatistic   deprecated
%   cfg.voxelthreshold   deprecated

% Copyright (C) 2005-2007, Robert Oostenveld
%
% $Log: statistics_montecarlo.m,v $
% Revision 1.22  2008/01/15 11:32:41  roboos
% implemented cfg.correctp, default is the same as before. Probably the
% old behaviour is UNDESIRED, but more discussion is needed and a thourough
% explanation should be provided on the mailing list and in the MEG meeting.
%
% Revision 1.21  2007/07/30 21:51:22  erimar
% Corrected a bug with respect to the copying of cfg.cvar in tmpcfg.cvar.
% In fact, cfg.cvar was not copied in tmpcfg.cvar. This is now corrected.
%
% Revision 1.20  2007/07/23 14:07:19  roboos
% fixed number of output arguments for statfun
%
% Revision 1.19  2007/07/23 13:26:52  jansch
% replaced feval by function-handle call
%
% Revision 1.18  2007/07/18 15:09:46  roboos
% Adjusted to the modified output of resampledesign, which now also returns
% the reindexing matrix for permutation (just like in case of bootstrap)
% instead of the permuted design itself.
% Renamed the "res" variable into "resample".
% Updated documentation with cfg.cvar.
%
% Revision 1.17  2007/07/17 10:37:06  roboos
% fixed bug for resampling dat in case of bootstrap (incorrect indexing)
% remember the cfg that is returned by clusterstat
%
% Revision 1.16  2007/07/04 16:04:53  roboos
% renamed randomizedesign to resampledesign
% first implementation of bootstrapping, sofar only for unpaired data
% added cfg.resampling=bootstrap/permutation
%
% Revision 1.15  2007/06/19 14:11:16  roboos
% only do backward compatibility for cfg.clusterthreshold if present
%
% Revision 1.14  2007/06/19 12:52:37  roboos
% implemented seperate common and individual nonparametric thresholds
% renamed the option cfg.clusterthreshold=nonparametric into nonparametric_individual, the new option is nonparametric_common
% updated documentation
%
% Revision 1.13  2007/03/27 15:38:00  erimar
% Passed cfg.dimord and cfg.dim to the called function. Updated help.
%
% Revision 1.12  2007/03/21 20:00:43  roboos
% Fixed bug (relates to cluster thresholding): also pass the clustertail
% to the statfun for computing the parametric clustercritval. The incorrect
% default was to use clustertail=0, also in case of one-tailed tests.
%
% Revision 1.11  2006/11/23 11:08:19  roboos
% implemented Holms' correction method for multiple comparisons
%
% Revision 1.10  2006/10/05 10:20:52  roboos
% updated cdocumentation
%
% Revision 1.9  2006/09/12 12:13:57  roboos
% fixed bug related to computation of clustercritval, ivar and uvar are needed for that
%
% Revision 1.8  2006/07/27 07:58:27  roboos
% improved documentation
%
% Revision 1.7  2006/07/20 15:35:47  roboos
% fixed a bug in the handling of some old cfgs
% added support for multi-field stat-structures as output of the statfun (only for statobs)
%
% Revision 1.6  2006/07/14 11:03:18  roboos
% remove clustertail from cfg in case no clustering is done
%
% Revision 1.5  2006/07/14 06:36:37  roboos
% changed handling of old cfg.correctm specification
%
% Revision 1.4  2006/07/14 06:29:29  roboos
% return cfg as output for later reference
%
% Revision 1.3  2006/07/12 15:53:25  roboos
% Moved all code and documentation from statistics_random to statistics_montecarlo as agreed with Eric
% The old functino is only there for backward compatibility, it gives an instructive warning
%
% Revision 1.2  2006/06/14 11:57:54  roboos
% changed documentation for cfg.correctm
%
% Revision 1.1  2006/05/31 13:04:31  roboos
% new implementation
%

% set the defaults for the main function
if ~isfield(cfg, 'alpha'),               cfg.alpha = 0.05;               end
if ~isfield(cfg, 'tail'),                cfg.tail = 0;                   end
if ~isfield(cfg, 'correctm'),            cfg.correctm = 'no';            end % no, max, cluster, bonferoni, fdr
if ~isfield(cfg, 'resampling'),          cfg.resampling = 'permutation'; end % permutation, bootstrap
if ~isfield(cfg, 'feedback'),            cfg.feedback = 'text';          end
if ~isfield(cfg, 'ivar'),                cfg.ivar = 'all';               end
if ~isfield(cfg, 'uvar'),                cfg.uvar = [];                  end
if ~isfield(cfg, 'cvar'),                cfg.cvar = [];                  end
if ~isfield(cfg, 'wvar'),                cfg.wvar = [];                  end
if ~isfield(cfg, 'correctp'),            cfg.correctp = 'no';            end % for the number of tails in a two-sided test

% for backward compatibility with old cfgs
if isfield(cfg, 'clusterstatistic') && ~isempty(cfg.clusterstatistic) && ~strcmp(cfg.clusterstatistic, 'no')
  warning('using cluster-based statistic for multiple comparison correction');
  cfg.correctm = 'cluster';
elseif strcmp(cfg.correctm, 'yes')
  warning('using maximum statistic for multiple comparison correction');
  cfg.correctm = 'max';
end

% set the defaults for clustering
if strcmp(cfg.correctm, 'cluster')
  if ~isfield(cfg, 'clusterstatistic'),    cfg.clusterstatistic = 'maxsum';     end  % no, max, maxsize, maxsum, wcm
  if ~isfield(cfg, 'clusterthreshold'),    cfg.clusterthreshold = 'parametric'; end  % parametric, nonparametric_individual, nonparametric_common
  if ~isfield(cfg, 'clusteralpha'),        cfg.clusteralpha = 0.05;             end
  if ~isfield(cfg, 'clustercritval'),      cfg.clustercritval = [];             end
  if ~isfield(cfg, 'clustertail'),         cfg.clustertail = cfg.tail;          end
else
  try, cfg = rmfield(cfg, 'clusterstatistic'); end
  try, cfg = rmfield(cfg, 'clusteralpha');     end
  try, cfg = rmfield(cfg, 'clustercritval');   end
  try, cfg = rmfield(cfg, 'clusterthreshold'); end
  try, cfg = rmfield(cfg, 'clustertail');      end
end

% for backward compatibility
if isfield(cfg, 'clusterthreshold') && strcmp(cfg.clusterthreshold, 'nonparametric')
  cfg.clusterthreshold = 'nonparametric_individual';
end

if size(design,2)~=size(dat,2)
  design = transpose(design);
end

if isfield(cfg, 'factor')
  cfg.ivar = cfg.factor;
  cfg = rmfield(cfg, 'factor');
end

if isfield(cfg, 'unitfactor')
  cfg.uvar = cfg.unitfactor;
  cfg = rmfield(cfg, 'unitfactor');
elseif isfield(cfg, 'repeatedmeasures')
  cfg.uvar = cfg.repeatedmeasures;
  cfg = rmfield(cfg, 'repeatedmeasures');
end

if ischar(cfg.ivar) && strcmp(cfg.ivar, 'all')
  cfg.ivar = 1:size(design,1);
end

% for backward compatibility in some older statfuns
cfg.factor     = cfg.ivar;
cfg.unitfactor = cfg.uvar;

% for backward compatibility with some low-level statfuns
if strcmp(cfg.statistic, 'corrcoef'),      cfg.statistic='pearson';      end
if strcmp(cfg.statistic, 'difference'),    cfg.statistic='diff';         end
if strcmp(cfg.statistic, 'anova'),         cfg.statistic='fstat';        end
if strcmp(cfg.statistic, 'paired-tstat'),  cfg.statistic='paired_tstat'; end

% for backward compatibility in clustering
if isfield(cfg, 'voxelthreshold') && strcmp(cfg.voxelstatistic, 'prob')
  cfg.clusteralpha     = cfg.voxelthreshold;
  cfg.clusterthreshold = 'nonparametric_individual';
  cfg.clustercritval   = [];
  cfg = rmfield(cfg, 'voxelthreshold');
  cfg = rmfield(cfg, 'voxelstatistic');
elseif isfield(cfg, 'voxelthreshold') && ~strcmp(cfg.voxelstatistic, 'prob')
  cfg.clusteralpha     = [];
  cfg.clusterthreshold = 'parametric';
  cfg.clustercritval   = cfg.voxelthreshold;
  cfg = rmfield(cfg, 'voxelthreshold');
  cfg = rmfield(cfg, 'voxelstatistic');
end

% the backward compatibility is broken here
if isfield(cfg, 'ztransform')
  error('cfg.ztransform is not supported any more');
end
if isfield(cfg, 'removemarginalmeans')
  error('cfg.removemarginalmeans is not supported any more');
end
if isfield(cfg, 'randomfactor'),
  error('cfg.randomfactor is not supported any more');
end

% determine the function handle to the low-level statistics function
if exist(['statistics_' cfg.statistic])
  statfun = str2func(['statistics_' cfg.statistic]);
elseif exist(['statfun_' cfg.statistic])
  statfun = str2func(['statfun_' cfg.statistic]);
else
  error(sprintf('could not find the statistics function "%s"\n', ['statfun_' cfg.statistic]));
end
fprintf('using "%s" for the single-sample statistics\n', func2str(statfun));

% construct the resampled design matrix or data-shuffling matrix
fprintf('constructing randomized design\n');
resample = resampledesign(cfg, design);
Nrand = size(resample,1);

% most of the statfuns result in this warning, which is not interesting
warning('off', 'MATLAB:warn_r14_stucture_assignment');

if strcmp(cfg.correctm, 'cluster')
  % determine the critical value for cluster thresholding
  if strcmp(cfg.clusterthreshold, 'nonparametric_individual') || strcmp(cfg.clusterthreshold, 'nonparametric_common')
    fprintf('using a nonparmetric threshold for clustering\n');
    cfg.clustercritval = [];  % this will be determined later
  elseif strcmp(cfg.clusterthreshold, 'parametric') && isempty(cfg.clustercritval)
    fprintf('computing a parmetric threshold for clustering\n');
    tmpcfg = [];
    tmpcfg.dimord         = cfg.dimord;
    tmpcfg.dim            = cfg.dim;
    tmpcfg.alpha          = cfg.clusteralpha;
    tmpcfg.tail           = cfg.clustertail;
    tmpcfg.design         = cfg.design;
    tmpcfg.ivar           = cfg.ivar;
    tmpcfg.uvar           = cfg.uvar;
    tmpcfg.cvar           = cfg.cvar;
    tmpcfg.wvar           = cfg.wvar;
    tmpcfg.computecritval = 'yes';  % explicitly request the computation of the crtitical value
    tmpcfg.computestat    = 'no';   % skip the computation of the statistic
    try
      cfg.clustercritval    = getfield(statfun(tmpcfg, dat, design), 'critval');
    catch
      disp(lasterr);
      error('could not determine the parametric critical value for clustering');
    end
  elseif strcmp(cfg.clusterthreshold, 'parametric') && ~isempty(cfg.clustercritval)
    fprintf('using the specified parametric threshold for clustering\n');
    cfg.clusteralpha = [];
  end
end

% compute the statistic for the observed data
progress('init', cfg.feedback, 'computing statistic');
% get an estimate of the time required per evaluation of the statfun
time_pre = cputime;

try
  % the nargout function in Matlab 6.5 and older does not work on function handles
  num = nargout(statfun);
catch
  num = 1;
end
if num>1
  % both the statistic and the (updated) configuration are returned
  [statobs, cfg] = statfun(cfg, dat, design);
else
  % only the statistic is returned
  statobs = statfun(cfg, dat, design);
end

if isstruct(statobs)
  % remember all details for later reference, continue to work with the statistic
  statfull = statobs;
  statobs  = getfield(statfull, 'stat');
else
  % remember the statistic for later reference, continue to work with the statistic
  statfull.stat = statobs;
end
time_eval = cputime - time_pre;
fprintf('estimated time per randomization is %d seconds\n', round(time_eval));

% pre-allocate some memory
if strcmp(cfg.correctm, 'cluster')
  statrand = zeros(size(statobs,1), size(resample,1));
else
  prb_pos   = zeros(size(statobs));
  prb_neg   = zeros(size(statobs));
end

% compute the statistic for the randomized data and count the outliers
for i=1:Nrand
  progress(i/Nrand, 'computing statistic %d from %d\n', i, Nrand);
  if strcmp(cfg.resampling, 'permutation')
    tmpdesign = design(:,resample(i,:));     % the columns in the design matrix are reshufled by means of permutation
    tmpdat    = dat;                        % the data itself is not shuffled
    if size(tmpdesign,1)==size(tmpdat,2)
      tmpdesign = transpose(tmpdesign);
    end
  elseif strcmp(cfg.resampling, 'bootstrap')
    tmpdesign = design;                     % the design matrix is not shuffled
    tmpdat    = dat(:,resample(i,:));        % the columns of the data are resampled by means of bootstrapping
  end
  if strcmp(cfg.correctm, 'cluster')
    % keep each randomization in memory for cluster postprocessing
    dum = statfun(cfg, tmpdat, tmpdesign);
    if isstruct(dum)
      statrand(:,i) = getfield(dum, 'stat');
    else
      statrand(:,i) = dum;
    end
  else
    % do not keep each randomization in memory, but process them on the fly
    statrand = statfun(cfg, tmpdat, tmpdesign);
    if isstruct(statrand)
      statrand = getfield(statrand, 'stat');
    end
    % the following line is for debugging 
    % stat.statkeep(:,i) = statrand;
    if strcmp(cfg.correctm, 'max')
      % compare each data element with the maximum statistic
      prb_pos = prb_pos + (statobs<max(statrand(:)));
      prb_neg = prb_neg + (statobs>min(statrand(:)));
      ref_pos(i) = max(statrand(:));
      ref_neg(i) = min(statrand(:));
    else
      % compare each data element with its own statistic
      prb_pos = prb_pos + (statobs<statrand);
      prb_neg = prb_neg + (statobs>statrand);
    end
  end
end
progress('close');

if strcmp(cfg.correctm, 'cluster')
  % do the cluster postprocessing
  [stat, cfg] = clusterstat(cfg, statrand, statobs);
else
  switch cfg.tail
  case 1
    clear prb_neg  % not needed any more, free some memory
    stat.prob = prb_pos./Nrand;
  case -1
    clear prb_pos  % not needed any more, free some memory
    stat.prob = prb_neg./Nrand;
  case 0
    % for each observation select the tail that corresponds with the lowest probability
    prb_neg = prb_neg./Nrand;
    prb_pos = prb_pos./Nrand;
    stat.prob = min(prb_neg, prb_pos); % this is the probability for the most unlikely tail
  end
end

% In case of a two tailed test, the type-I error rate (alpha) refers to
% both tails of the distribution, whereas the stat.prob value computed sofar
% corresponds with one tail, i.e. the probability, under the assumption of
% no effect or no difference (the null hypothesis), of obtaining a result
% equal to or more extreme than what was actually observed. The decision
% rule whether the null-hopothesis should be rejected given the observed
% probability therefore should consider alpha divided by two, to correspond
% with the probability in one of the tails (the most extreme tail). This
% is conceptually equivalent to performing a Bonferoni correction for the
% two tails.
% 
% An alternative solution to distributed the alpha level over both tails is
% achieved by multiplying the probability with a factor of two, prior to
% thresholding it wich cfg.alpha.  The advantage of this solution is that
% it results in a p-value that corresponds with a parametric probability.
if strcmp(cfg.correctp, 'yes') && cfg.tail==0
  stat.prob = stat.prob .* 2;
end

switch lower(cfg.correctm)
  case 'max'
    % the correction is implicit in the method
    fprintf('using a maximum-statistic based method for multiple comparison correction\n');
    fprintf('the returned probabilities and the thresholded mask are corrected for multiple comparisons\n');
    stat.mask = stat.prob<=cfg.alpha;
    stat.ref_pos = ref_pos;
    stat.ref_neg = ref_neg;
  case 'cluster'
    % the correction is implicit in the method
    fprintf('using a cluster-based method for multiple comparison correction\n');
    fprintf('the returned probabilities and the thresholded mask are corrected for multiple comparisons\n');
    stat.mask = stat.prob<=cfg.alpha;
  case 'bonferoni'
    fprintf('performing Bonferoni correction for multiple comparisons\n');
    fprintf('the returned probabilities are uncorrected, the thresholded mask is corrected\n');
    stat.mask = stat.prob<=(cfg.alpha ./ numel(stat.prob));
  case 'holms'
    % test the most significatt significance probability against alpha/N, the second largest against alpha/(N-1), etc.
    fprintf('performing Holms correction for multiple comparisons\n');
    fprintf('the returned probabilities are uncorrected, the thresholded mask is corrected\n');
    [pvals, indx] = sort(stat.prob(:));                     % this sorts the significance probabilities from smallest to largest
    mask = pvals<=(cfg.alpha ./ ((length(pvals):-1:1)'));   % compare each significance probability against its individual threshold
    stat.mask = zeros(size(stat.prob));
    stat.mask(indx) = mask;
  case 'fdr'
    fprintf('performing FDR correction for multiple comparisons\n');
    fprintf('the returned probabilities are uncorrected, the thresholded mask is corrected\n');
    stat.mask = fdr(stat.prob, cfg.alpha);
  otherwise
    fprintf('not performing a correction for multiple comparisons\n');
    stat.mask = stat.prob<=cfg.alpha;
end

% return the observed statistic
if ~isfield(stat, 'stat')
  stat.stat = statobs;
end

% return optional other details that were returned by the statfun
fn = fieldnames(statfull);
for i=1:length(fn)
  if ~isfield(stat, fn{i})
    stat = setfield(stat, fn{i}, getfield(statfull, fn{i}));
  end
end
