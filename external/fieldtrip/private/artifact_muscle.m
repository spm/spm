function [cfg, artifact] = artifact_muscle(cfg)

% ARTIFACT_MUSCLE reads the data segments of interest from file and
% identifies muscle artifacts.
%
% Use as
%   [cfg, artifact] = artifact_muscle(cfg)
% where the configuration structure should contain
%   cfg.trl = structure that defines the data segments of interest. See DEFINETRIAL
%   cfg.dataset    = string, name of the dataset
% or, instead of cfg.dataset
%   cfg.datafile   = string
%   cfg.headerfile = string
%
% The raw data is preprocessed with the following configuration parameters,
% which are optimal for identifying muscle artifacts:
%   cfg.artfctdef.muscle.bpfilter    = 'yes'
%   cfg.artfctdef.muscle.bpfreq      = [110 140]
%   cfg.artfctdef.muscle.bpfiltord   = 10
%   cfg.artfctdef.muscle.bpfilttype  = 'but'
%   cfg.artfctdef.muscle.hilbert     = 'yes'
%   cfg.artfctdef.muscle.boxcar      = 0.2
%
% Artifacts are identified by means of thresholding the z-transformed value
% of the preprocessed data.
%   cfg.artfctdef.muscle.channel     = Nx1 cell-array with selection of channels, see CHANNELSELECTION for details
%   cfg.artfctdef.muscle.cutoff      = 4       z-value at which to threshold
%   cfg.artfctdef.muscle.trlpadding  = 0.1
%   cfg.artfctdef.muscle.fltpadding  = 0.1
%   cfg.artfctdef.muscle.artpadding  = 0.1
%
% The output argument "artifact" is a Nx2 matrix comparable to the
% "trl" matrix of DEFINETRIAL. The first column of which specifying the
% beginsamples of an artifact period, the second column contains the
% endsamples of the artifactperiods.
%
% See also ARTIFACT_ZVALUE, REJECTARTIFACT
  
% Undocumented local options:
% cfg.datatype
% cfg.method
%
% This function depends on ARTIFACT_ZVALUE which has the following options:
% cfg.artfctdef.zvalue.channel in ARTIFACT_MUSCLE as cfg.artfctdef.muscle.channel default = 'MEG', documented
% cfg.artfctdef.zvalue.cutoff in ARTIFACT_MUSCLE as cfg.artfctdef.muscle.cutoff default = 4, documented
% cfg.artfctdef.zvalue.trlpadding in ARTIFACT_MUSCLE as cfg.artfctdef.muscle.trlpadding default = 0.1, documented
% cfg.artfctdef.zvalue.fltpadding in ARTIFACT_MUSCLE as cfg.artfctdef.muscle.fltpadding default = 0.1, documented
% cfg.artfctdef.zvalue.artpadding in ARTIFACT_MUSCLE as cfg.artfctdef.muscle.artpadding default = 0.1, documented
% cfg.artfctdef.zvalue.feedback


% Copyright (c) 2003-2006, Jan-Mathijs Schoffelen & Robert Oostenveld
%
% $Log: artifact_muscle.m,v $
% Revision 1.21  2006/12/04 09:30:02  roboos
% changed a comment
%
% Revision 1.20  2006/11/29 09:06:36  roboos
% renamed all cfg options with "sgn" into "channel", added backward compatibility when required
% updated documentation, mainly in the artifact detection routines
%
% Revision 1.19  2006/04/25 17:06:28  ingnie
% updated documentation
%
% Revision 1.18  2006/04/20 09:58:33  roboos
% updated documentation
%
% Revision 1.17  2006/02/27 17:01:07  roboos
% added tmpcfg.dataset, added try-end around copying of headerfile, datafile and dataset
%
% Revision 1.16  2006/01/12 13:51:38  roboos
% completely new implementation, all based upon the same artifact_zvalue code
% all preprocessing is now done consistently and the various paddings have been better defined
% the functions do not have any explicit support any more for non-continuous data
% the old artifact_xxx functions from JM have been renamed to xxx_old
%

% set default rejection parameters
if ~isfield(cfg,'artfctdef'),                     cfg.artfctdef                    = [];        end
if ~isfield(cfg.artfctdef,'muscle'),              cfg.artfctdef.muscle             = [];        end
if ~isfield(cfg.artfctdef.muscle,'method'),       cfg.artfctdef.muscle.method      = 'zvalue';  end

% for backward compatibility
if isfield(cfg.artfctdef.muscle,'sgn')
  cfg.artfctdef.muscle.channel = cfg.artfctdef.muscle.sgn;
  cfg.artfctdef.muscle         = rmfield(cfg.artfctdef.muscle, 'sgn');
end

if isfield(cfg.artfctdef.muscle, 'artifact')
  fprintf('muscle artifact detection has already been done, retaining artifacts\n');
  artifact = cfg.artfctdef.muscle.artifact;
  return
end

if strcmp(cfg.artfctdef.muscle.method, 'zvalue')
  % the following settings should be supported for backward compatibility
  if isfield(cfg.artfctdef.muscle,'pssbnd'),
    cfg.artfctdef.muscle.bpfreq   = cfg.artfctdef.muscle.pssbnd;
    cfg.artfctdef.muscle.bpfilter = 'yes';
    cfg.artfctdef.muscle = rmfield(cfg.artfctdef.muscle,'pssbnd');
  end;
  dum = 0;
  if isfield(cfg.artfctdef.muscle,'pretim'),
    dum = max(dum, cfg.artfctdef.muscle.pretim);
    cfg.artfctdef.muscle = rmfield(cfg.artfctdef.muscle,'pretim');
  end
  if isfield(cfg.artfctdef.muscle,'psttim'),
    dum = max(dum, cfg.artfctdef.muscle.psttim);
    cfg.artfctdef.muscle = rmfield(cfg.artfctdef.muscle,'psttim');
  end
  if dum
    cfg.artfctdef.muscle.artpadding = max(dum);
  end
  if isfield(cfg.artfctdef.muscle,'padding'),
    cfg.artfctdef.muscle.trlpadding   = cfg.artfctdef.muscle.padding;
    cfg.artfctdef.muscle = rmfield(cfg.artfctdef.muscle,'padding');
  end

  % settings for preprocessing
  if ~isfield(cfg.artfctdef.muscle,'bpfilter'),   cfg.artfctdef.muscle.bpfilter    = 'yes';     end
  if ~isfield(cfg.artfctdef.muscle,'bpfreq'),     cfg.artfctdef.muscle.bpfreq      = [110 140]; end
  if ~isfield(cfg.artfctdef.muscle,'bpfiltord'),  cfg.artfctdef.muscle.bpfiltord   = 10;        end
  if ~isfield(cfg.artfctdef.muscle,'bpfilttype'), cfg.artfctdef.muscle.bpfilttype  = 'but';     end
  if ~isfield(cfg.artfctdef.muscle,'hilbert'),    cfg.artfctdef.muscle.hilbert     = 'yes';     end
  if ~isfield(cfg.artfctdef.muscle,'boxcar'),     cfg.artfctdef.muscle.boxcar      = 0.2;       end
  % settings for the zvalue subfunction
  if ~isfield(cfg.artfctdef.muscle,'channel'),    cfg.artfctdef.muscle.channel     = 'MEG';     end
  if ~isfield(cfg.artfctdef.muscle,'trlpadding'), cfg.artfctdef.muscle.trlpadding  = 0.1;       end
  if ~isfield(cfg.artfctdef.muscle,'fltpadding'), cfg.artfctdef.muscle.fltpadding  = 0.1;       end
  if ~isfield(cfg.artfctdef.muscle,'artpadding'), cfg.artfctdef.muscle.artpadding  = 0.1;       end
  if ~isfield(cfg.artfctdef.muscle,'cutoff'),     cfg.artfctdef.muscle.cutoff      = 4;         end
  % construct a temporary configuration that can be passed onto artifact_zvalue
  tmpcfg                  = [];
  tmpcfg.trl              = cfg.trl;
  try, tmpcfg.dataset     = cfg.dataset; end
  try, tmpcfg.datafile    = cfg.datafile; end
  try, tmpcfg.headerfile  = cfg.headerfile; end
  tmpcfg.artfctdef.zvalue = cfg.artfctdef.muscle;
  if isfield(cfg, 'datatype')
    tmpcfg.datatype = cfg.datatype;
  end
  % call the zvalue artifact detection function
  [tmpcfg, artifact]  = artifact_zvalue(tmpcfg);
  cfg.artfctdef.muscle = tmpcfg.artfctdef.zvalue;
else
  error(sprintf('muscle artifact detection only works with cfg.method=''zvalue'''));
end

