function [comp] = componentanalysis(cfg, data);

% COMPONENTANALYSIS principal or independent component analysis
% computes the topography and timecourses of the ICA/PCA components 
% in the EEG/MEG data.
%
% Use as
%   [comp] = componentanalysis(cfg, data)
%
% where the data comes from PREPROCESING or TIMELOCKANALYSIS and the
% configuration structure can contain
%   cfg.method       = 'runica', 'binica', 'pca', 'jader', 'varimax', 'dss', 'cca' (default = 'runica')
%   cfg.channel      = cell-array with channel selection (default = 'all'), see CHANNELSELECTION for details
%   cfg.trials       = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.numcomponent = 'all' or number (default = 'all')
%   cfg.blc          = 'no' or 'yes' (default = 'yes')
%   cfg.detrend      = 'no' or 'yes' (default = 'no')
%   cfg.runica       = substructure with additional low-level options for this method
%   cfg.binica       = substructure with additional low-level options for this method
%   cfg.dss          = substructure with additional low-level options for this method
%   cfg.fastica      = substructure with additional low-level options for this method
% 
% Instead of specifying a component analysis method, you can also specify
% a previously computed mixing matrix, which will be used to estimate the
% component timecourses in this data. This requires
%   cfg.topo         = NxN matrix with a component topography in each column
%   cfg.topolabel    = Nx1 cell-array with the channel labels
% 
% See also FASTICA, RUNICA, SVD, JADER, VARIMAX, DSS, CCA

% NOTE parafac is also implemented, but that does not fit into the
% structure of 2D decompositions very well. Probably I should implement it
% in a separate function for N-D decompositions

% Copyright (C) 2003-2007, Robert Oostenveld
%
% $Log: componentanalysis.m,v $
% Revision 1.34  2008/06/17 15:21:42  sashae
% now using preproc_modules
%
% Revision 1.33  2008/05/06 14:23:32  sashae
% change in trial selection, cfg.trials can be a logical
%
% Revision 1.32  2007/12/18 17:06:05  sashae
% updated documentation
%
% Revision 1.31  2007/05/02 15:59:13  roboos
% be more strict on the input and output data: It is now the task of
% the private/checkdata function to convert the input data to raw
% data (i.e. as if it were coming straight from preprocessing).
% Furthermore, the output data is NOT converted back any more to the
% input data, i.e. the output data is the same as what it would be
% on raw data as input, regardless of the actual input.
%
% Revision 1.30  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.29  2007/03/30 17:05:40  ingnie
% checkinput; only proceed when input data is allowed datatype
%
% Revision 1.28  2007/03/04 14:22:14  chrhes
% changes to pca option: this is now calculated using eigenvalue decomposition
% (EVD) of the data cross-covariance matrix (the traditional approach to PCA)
% rather than by means of singular value decomposition (SVD) of the data matrix
% (which has been added as a separate option: svd). For large data matrices
% EVD of the cross-covariancre matrix is computationally much faster than SVD,
% and avoids the risk of memory related errors in matlab.
%
% Revision 1.27  2007/03/04 13:45:12  chrhes
% changes to how the code calls fastica: the responsibility for providing the 
% correct optional arguments through cfg.fastica lies entriely with the user 
% (as for the runica option)
%
% Revision 1.26  2007/02/27 09:54:22  roboos
% added required defaults for binica, check for EEGLAB in case of binica
%
% Revision 1.25  2007/02/12 19:57:28  roboos
% implemented binica
%
% Revision 1.24  2007/02/12 19:44:31  roboos
% check fastica toolbox, try to add automatically
%
% Revision 1.23  2006/12/15 14:44:40  chrhes
% implemented (undocumented, preliminary) support for the FastICA algorithm
%
% Revision 1.22  2006/10/04 07:10:07  roboos
% updated documentation
%
% Revision 1.21  2006/08/16 10:49:53  roboos
% updated documentation
%
% Revision 1.20  2006/08/16 10:46:53  roboos
% added support for unmixing data using a previously determined (un)mixing 
% matrix (options cfg.topo and cfg.topolabel) channels such as ECG can be 
% excluded from teh unmixing, but will remain in the output data
%
% Revision 1.19  2006/06/07 09:34:19  roboos
% changed checktoolbox into hastoolbox
%
% Revision 1.18  2006/04/20 09:58:33  roboos
% updated documentation
%
% Revision 1.17  2006/01/11 12:58:16  roboos
% minor change to documentation
%
% Revision 1.16  2006/01/06 11:37:13  roboos
% switched to the checktoolbox() subfunction for toolbox dependency checking
% implemented canonical correlation analysis (CCA) using a function from Magnus 
% Borga changed the DSS implementation from version 0.6 beta to version 1.0
% fixed some small bugs related to the number of components
%
% Revision 1.15  2005/11/08 11:36:35  roboos
% updated the documentation to include dss
%
% Revision 1.14  2005/11/04 17:14:41  roboos
% implemented support for denoising source separation (dss), requires external 
% toolbox
%
% Revision 1.13  2005/11/04 10:26:03  roboos
% changed the handling of optional arguments for runica to use cfg2keyval()
%
% Revision 1.12  2005/10/18 12:23:46  roboos
% added tic/toc timer around complete function
% changed handing of averages to using data2raw and raw2data
% added support for more options of runica, using cfg.runica substructure
%
% Revision 1.11  2005/08/05 09:16:22  roboos
% removed the obsolete data.offset
%
% Revision 1.10  2005/06/29 12:46:29  roboos
% the following changes were done on a large number of files at the same time, 
% not all of them may apply to this specific file
% - added try-catch around the inclusion of input data configuration
% - moved cfg.version, cfg.previous and the assignment of the output cfg to the end
% - changed help comments around the configuration handling
% - some changes in whitespace
%
% Revision 1.9  2005/06/28 16:08:00  roboos
% cleaned up the detectin and error handling for eeglab and nway toolbox dependencies
%
% Revision 1.8  2004/09/16 15:36:08  roboos
% added preliminary support for parafac
%
% Revision 1.7  2004/09/01 17:59:28  roboos
% added copyright statements to all filed
% added cfg.version to all functions that give configuration in their output
% added cfg.previous to all functions with input data containing configuration details
%
% Revision 1.6  2004/05/19 14:47:20  roberto
% added version details to output configuration
%
% Revision 1.5  2004/03/06 13:06:49  roberto
% cleaded up transposing of weights matrix, no functional change
%
% Revision 1.4  2004/02/05 10:31:28  roberto
% implemented support for erf/erp input (from timelockanalysis)
% added option for svds with reduced number of components
% some layout changes to the code and updated documentation
%
% Revision 1.3  2004/02/04 14:01:44  roberto
% only change in help
%
% Revision 1.2  2003/12/16 21:03:19  roberto
% renamed my_xxx to their original function names (my_runica->runica etc.)
%
% Revision 1.1  2003/11/12 07:57:23  roberto
% new implementation, works with two functions from eeglab
%

% set a timer to determine how long this function takes
tic;

% check if the input data is valid for this function
data = checkdata(data, 'datatype', 'raw', 'feedback', 'yes');

% set the defaults
if ~isfield(cfg, 'method'),        cfg.method  = 'runica';     end
if ~isfield(cfg, 'blc'),           cfg.blc     = 'yes';        end
if ~isfield(cfg, 'detrend'),       cfg.detrend = 'no';         end
if ~isfield(cfg, 'trials'),        cfg.trials  = 'all';        end
if ~isfield(cfg, 'channel'),       cfg.channel = 'all';        end
if ~isfield(cfg, 'numcomponent'),  cfg.numcomponent = 'all';   end

if isfield(cfg, 'topo') && isfield(cfg, 'topolabel')
  % use the previously determined unmixing matrix on this dataset
  % remove all cfg settings  that do not apply 
  tmpcfg = [];
  tmpcfg.blc       = cfg.blc;
  tmpcfg.detrend   = cfg.detrend;
  tmpcfg.trials    = cfg.trials;
  tmpcfg.channel   = cfg.channel;
  tmpcfg.topo      = cfg.topo;
  tmpcfg.topolabel = cfg.topolabel;
  tmpcfg.numcomponent = 'all';
  tmpcfg.method    = 'predetermined mixing matrix';
  cfg = tmpcfg;
end

% additional options, see FASTICA for details
if ~isfield(cfg, 'fastica'),        cfg.fastica = [];          end;

% additional options, see RUNICA for details
if ~isfield(cfg, 'runica'),        cfg.runica = [];            end
if ~isfield(cfg.runica, 'lrate'),  cfg.runica.lrate = 0.001;   end

% additional options, see BINICA for details
if ~isfield(cfg, 'binica'),        cfg.binica = [];            end
if ~isfield(cfg.binica, 'lrate'),  cfg.binica.lrate = 0.001;   end

% additional options, see DSS for details
if ~isfield(cfg, 'dss'),                  cfg.dss      = [];                           end
if ~isfield(cfg.dss, 'denf'),             cfg.dss.denf = [];                           end
if ~isfield(cfg.dss.denf, 'function'),    cfg.dss.denf.function = 'denoise_fica_tanh'; end
if ~isfield(cfg.dss.denf, 'params'),      cfg.dss.denf.params   = [];                  end

% check whether the required low-level toolboxes are installed
switch cfg.method
   case 'fastica'
      hastoolbox('fastica', 1);       % see http://www.cis.hut.fi/projects/ica/fastica
   case {'runica', 'jader', 'varimax', 'binica'}
      hastoolbox('eeglab', 1);        % see http://www.sccn.ucsd.edu/eeglab
   case 'parafac'
      hastoolbox('nway', 1);          % see http://www.models.kvl.dk/source/nwaytoolbox
   case 'dss'
      hastoolbox('dss', 1);           % see http://www.cis.hut.fi/projects/dss
end % cfg.method

% default is to compute just as many components as there are channels in the data
if strcmp(cfg.numcomponent, 'all')
  cfg.numcomponent = length(data.label);
end

% select trials of interest
if ~strcmp(cfg.trials, 'all')
  if islogical(cfg.trials),  cfg.trials=find(cfg.trials);  end
  fprintf('selecting %d trials\n', length(cfg.trials));
  data.trial  = data.trial(cfg.trials);
  data.time   = data.time(cfg.trials);
end
Ntrials  = length(data.trial);

% select channels of interest
cfg.channel = channelselection(cfg.channel, data.label);
chansel = match_str(data.label, cfg.channel);
fprintf('selecting %d channels\n', length(chansel));
for trial=1:Ntrials
  data.trial{trial} = data.trial{trial}(chansel,:);
end
data.label = data.label(chansel);
Nchans     = length(chansel);

% determine the size of each trial, they can be variable length
for trial=1:Ntrials
  Nsamples(trial) = size(data.trial{trial},2);
end

if strcmp(cfg.detrend, 'yes')
  % optionally perform detrending on each trial
  fprintf('detrending data \n');
  for trial=1:Ntrials
    data.trial{trial} = preproc_detrend(data.trial{trial});
  end
elseif strcmp(cfg.blc, 'yes') 
  % optionally perform baseline correction on each trial
  fprintf('baseline correcting data \n');
  for trial=1:Ntrials
    data.trial{trial} = preproc_baselinecorrect(data.trial{trial});
  end
end

if strcmp(cfg.method, 'predetermined mixing matrix')
  % the single trial data does not have to be concatenated
elseif strcmp(cfg.method, 'parafac')
  % concatenate all the data into a 3D matrix
  fprintf('concatenating data');
  Nsamples = Nsamples(1);
  dat = zeros(Ntrials, Nchans, Nsamples);
  % all trials should have an equal number of samples
  % and it is assumed that the time axes of all trials are aligned
  for trial=1:Ntrials
    fprintf('.');
    dat(trial,:,:) = data.trial{trial};
  end
  fprintf('\n');
  fprintf('concatenated data matrix size %dx%dx%d\n', size(dat,1), size(dat,2), size(dat,3));
else
  % concatenate all the data into a 2D matrix
  fprintf('concatenating data');
  dat = zeros(Nchans, sum(Nsamples));
  for trial=1:Ntrials
    fprintf('.');
    begsample = sum(Nsamples(1:(trial-1))) + 1;
    endsample = sum(Nsamples(1:trial));
    dat(:,begsample:endsample) = data.trial{trial};
  end
  fprintf('\n');
  fprintf('concatenated data matrix size %dx%d\n', size(dat,1), size(dat,2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform the component analysis
fprintf('starting decomposition using %s\n', cfg.method);
switch cfg.method

   case 'fastica'
      % construct key-value pairs for the optional arguments
      optarg = cfg2keyval(cfg.fastica);
      [A, W] = fastica(dat, optarg{:});
      weights = W;
      sphere = eye(size(W,2));

   case 'runica'
      % construct key-value pairs for the optional arguments
      optarg = cfg2keyval(cfg.runica);
      [weights, sphere] = runica(dat, optarg{:});

   case 'binica'
      % construct key-value pairs for the optional arguments
      optarg = cfg2keyval(cfg.binica);
      [weights, sphere] = binica(dat, optarg{:});

   case 'jader'
      weights = jader(dat);
      sphere  = eye(size(weights, 2));

   case 'varimax'
      weights = varimax(dat);
      sphere  = eye(size(weights, 2));

   case 'cca'
      [y, w] = ccabss(dat);
      weights = w';
      sphere  = eye(size(weights, 2));

   case 'pca'
      % compute data cross-covariance matrix
      C = (dat*dat')./(size(dat,2)-1);
      % eigenvalue decomposition (EVD)
      [E,D] = eig(C);
      % sort eigenvectors in descending order of eigenvalues 
      d = cat(2,[1:1:Nchans]',diag(D));
      d = sortrows(d,[-2]);
      % return the desired number of principal components
      weights = E(:,d(1:cfg.numcomponent,1))';
      sphere = eye(size(weights,2));
      clear C D E d
    
   case 'svd'
      if cfg.numcomponent<Nchans
         % compute only the first components
         [u, s, v] = svds(dat, cfg.numcomponent);
         u(Nchans, Nchans) = 0;
      else
         % compute all components
         [u, s, v] = svd(dat, 0);
      end
      weights = u';
      sphere  = eye(size(weights, 2));

   case 'parafac'
      f = parafac(dat, cfg.numcomponent);

   case 'dss'
      params         = cfg.dss;
      params.denf.h  = str2func(cfg.dss.denf.function);
      if ~ischar(cfg.numcomponent)
         params.sdim = cfg.numcomponent;
      end
      % create the state
      state   = dss_create_state(dat, params);
      % increase the amount of information that is displayed on screen
      state.verbose = 3;
      % start the decomposition
      % state   = dss(state);  % this is for the DSS toolbox version 0.6 beta
      state   = denss(state);  % this is for the DSS toolbox version 1.0
      weights = state.W;
      sphere  = state.V;
      % remember the updated configuration details
      cfg.dss.denf      = state.denf;
      cfg.numcomponent  = state.sdim;

   case 'predetermined mixing matrix'
      [compsel, datsel] = match_str(cfg.topolabel, data.label);
      cfg.numcomponent = length(compsel);
      Nchan   = length(data.label);
      notsel  = setdiff(1:Nchan, datsel);
      sphere  = eye(Nchan,Nchan);
      weights = eye(Nchan,Nchan);
      weights(datsel,datsel) = inv(cfg.topo(compsel,compsel));

   otherwise
      error('unknown method for component analysis');
end % switch method


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collect the results
try, comp.fsample = data.fsample; end
try, comp.time    = data.time;    end

% compute the activations in each trial
for trial=1:Ntrials
  if strcmp(cfg.method, 'parafac')
    % FIXME, this is not properly supported yet
    comp.trial{trial} = [];
  else
    comp.trial{trial} = weights * sphere * data.trial{trial};
  end
end

if strcmp(cfg.method, 'predetermined mixing matrix')
  % the output may consist of both unmixed components and of original channels
  comp.label = data.label;
  for i=1:cfg.numcomponent
    comp.label{datsel(i)} = sprintf('component%03d', i);
  end
else
  % the output only contains unmixed components
  for chan=1:cfg.numcomponent
    % create a new label for each component
    comp.label{chan} = sprintf('%s%03d', cfg.method, chan);
  end
  comp.label     = comp.label(:);          % convert labels into column vector
end
comp.topolabel   = data.label(:);          % these labels describe the topographic components

if strcmp(cfg.method, 'parafac')
  comp.topo = f{2};                        % topography of each component
  comp.f1 = f{1};                          % FIXME, this is not properly supported yet
  comp.f2 = f{2};                          % FIXME, this is not properly supported yet
  comp.f3 = f{3};                          % FIXME, this is not properly supported yet
else
   if (size(weights,1)==size(weights,2))
      comp.topo = inv(weights * sphere);      % topography of each component
   else
      comp.topo = pinv(weights * sphere);     % fix to allow fewer sources than sensors
   end
end

% add the version details of this function call to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id   = '$Id: componentanalysis.m,v 1.34 2008/06/17 15:21:42 sashae Exp $';
% remember the configuration details of the input data
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output 
comp.cfg = cfg;

fprintf('total time in componentanalysis %.1f seconds\n', toc);

