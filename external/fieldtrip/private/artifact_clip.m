function [cfg, artifact] = artifact_clip(cfg);

% ARTIFACT_CLIP scans the data segments of interest for channels that
% clip. A clipping artifact is detected by the signal being completely
% flat for some time.
%
% Use as
%   [cfg, artifact] = artifact_clip(cfg)
% where the configuration should contain
%   cfg.artfctdef.clip.channel  = Nx1 cell-array with selection of channels, see CHANNELSELECTION for details
%   cfg.artfctdef.clip.pretim   = 0.000;  pre-artifact rejection-interval in seconds
%   cfg.artfctdef.clip.psttim   = 0.000;  post-artifact rejection-interval in seconds
%   cfg.artfctdef.clip.thresh   = 0.010;  minimum duration in seconds of a datasegment with consecutive identical samples to be considered as 'clipped'
%
% See also REJECTARTIFACT

% Undocumented local options:
% cfg.trl
% cfg.version
%
% This function depends on PREPROC which has the following options:
% cfg.absdiff
% cfg.blc
% cfg.blcwindow
% cfg.boxcar
% cfg.bpfilter
% cfg.bpfiltord
% cfg.bpfilttype
% cfg.bpfreq
% cfg.derivative
% cfg.detrend
% cfg.dftfilter
% cfg.dftfreq
% cfg.hilbert
% cfg.hpfilter
% cfg.hpfiltord
% cfg.hpfilttype
% cfg.hpfreq
% cfg.implicitref
% cfg.lnfilter
% cfg.lnfiltord
% cfg.lnfreq
% cfg.lpfilter
% cfg.lpfiltord
% cfg.lpfilttype
% cfg.lpfreq
% cfg.medianfilter
% cfg.medianfiltord
% cfg.rectify
% cfg.refchannel
% cfg.reref

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: artifact_clip.m,v $
% Revision 1.9  2006/11/29 09:06:36  roboos
% renamed all cfg options with "sgn" into "channel", added backward compatibility when required
% updated documentation, mainly in the artifact detection routines
%
% Revision 1.8  2006/06/14 12:43:48  roboos
% removed the documentation for cfg.lnfilttype, since that option is not supported by preproc
%
% Revision 1.7  2006/04/20 09:58:33  roboos
% updated documentation
%
% Revision 1.6  2006/04/12 08:38:01  ingnie
% updated documentation
%
% Revision 1.5  2006/01/31 13:49:29  jansch
% included dataset2files to ensure the presence of cfg.headerfile or cfg.datafile
% whenever needed
%
% Revision 1.4  2005/12/20 08:36:47  roboos
% add the artifact Nx2 matrix to the output configuration
% changed some indentation and white space, renamed a few variables
%
% Revision 1.3  2005/10/06 14:28:20  roboos
% added check for case when no artifacts are found
%
% Revision 1.2  2005/09/05 10:37:03  roboos
% added pretim and psttim option for extending the detected artifacts
% added copyright and log
%

% set default rejection parameters for clip artifacts if necessary.
if ~isfield(cfg,'artfctdef'),               cfg.artfctdef               = [];              end;
if ~isfield(cfg.artfctdef,'clip'),          cfg.artfctdef.clip          = [];              end;
if ~isfield(cfg.artfctdef.clip,'channel'),  cfg.artfctdef.clip.channel  = 'all';           end;
if ~isfield(cfg.artfctdef.clip,'thresh'),   cfg.artfctdef.clip.thresh   = 0.010;           end;
if ~isfield(cfg.artfctdef.clip,'pretim'),   cfg.artfctdef.clip.pretim   = 0.000;           end;
if ~isfield(cfg.artfctdef.clip,'psttim'),   cfg.artfctdef.clip.psttim   = 0.000;           end;

% for backward compatibility
if isfield(cfg.artfctdef.clip,'sgn')
  cfg.artfctdef.clip.channel = cfg.artfctdef.clip.sgn;
  cfg.artfctdef.clip         = rmfield(cfg.artfctdef.clip, 'sgn');
end

% start with an empty artifact list
artifact = [];

% read the header
cfg = dataset2files(cfg);
hdr = read_fcdc_header(cfg.headerfile);

% find the channel labels present in the data and their indices
label = channelselection(cfg.artfctdef.clip.channel, hdr.label);
sgnindx = match_str(hdr.label, label);

% make a local copy for convenience
artfctdef = cfg.artfctdef.clip;

ntrl = size(cfg.trl,1);
nsgn = length(sgnindx);
for trlop=1:ntrl
  fprintf('searching for clipping artifacts in trial %d\n', trlop);
  % read the data of this trial
  dat = read_fcdc_data(cfg.datafile, hdr, cfg.trl(trlop,1), cfg.trl(trlop,2), sgnindx);
  % apply filtering etc to the data
  datflt = preproc(dat, label, hdr.Fs, artfctdef, cfg.trl(trlop,3));
  % detect all samples that have the same value as the previous sample
  identical = (datflt(:,1:(end-1)) == datflt(:,2:end));
  % ensure that the number of samples does not change
  identical = [identical zeros(nsgn,1)];

  % determine the number of consecutively identical samples
  clip = zeros(size(dat));
  for sgnlop=1:length(sgnindx)
    up = find(diff([0 identical(sgnlop,:)], 1, 2)== 1);
    dw = find(diff([identical(sgnlop,:) 0], 1, 2)==-1);
    for k=1:length(up)
      clip(sgnlop,up(k):dw(k)) = dw(k)-up(k);
    end
  end
  % collapse over cannels
  clip = max(clip,[],1);

  % detect whether there are intervals in which the number of consecutive
  % identical samples is larger than the threshold
  thresh = (clip>=artfctdef.thresh*hdr.Fs);

  % remember the thresholded parts as artifacts
  artup = find(diff([0 thresh])== 1) + cfg.trl(trlop,1) - 1;
  artdw = find(diff([thresh 0])==-1) + cfg.trl(trlop,1) - 1;
  for k=1:length(artup)
    artifact(end+1,:) = [artup(k) artdw(k)];
  end
end

if ~isempty(artifact)
  % add the pretim and psttim to the detected artifacts
  artifact(:,1) = artifact(:,1) - artfctdef.pretim * hdr.Fs;
  artifact(:,2) = artifact(:,2) + artfctdef.psttim * hdr.Fs;
end

% remember the details that were used here
cfg.artfctdef.clip          = artfctdef;
cfg.artfctdef.clip.label    = label;
cfg.artfctdef.clip.trl      = cfg.trl;
cfg.artfctdef.clip.artifact = artifact;

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: artifact_clip.m,v 1.9 2006/11/29 09:06:36 roboos Exp $';
