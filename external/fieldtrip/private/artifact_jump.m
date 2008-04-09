function [cfg, artifact] = artifact_jump(cfg)

% ARTIFACT_JUMP reads the data segments of interest from file and
% identifies SQUID jump artifacts.
%
% Use as
%   [cfg, artifact] = artifact_jump(cfg)
% where the configuration structure should contain
%   cfg.trl = structure that defines the data segments of interest. See DEFINETRIAL
%   cfg.dataset    = string, name of the dataset
% or, instead of cfg.dataset
%   cfg.datafile   = string
%   cfg.headerfile = string
%
% The raw data is preprocessed with the following configuration parameters,
% which are optimal for identifying jump artifacts:
%   cfg.artfctdef.jump.medianfilter  = 'yes'
%   cfg.artfctdef.jump.medianfiltord = 9
%   cfg.artfctdef.jump.absdiff       = 'yes'
%
% Artifacts are identified by means of thresholding the z-transformed value
% of the preprocessed data.
%   cfg.artfctdef.jump.channel       = Nx1 cell-array with selection of channels, see CHANNELSELECTION for details
%   cfg.artfctdef.jump.cutoff        = 20      z-value at which to threshold
%   cfg.artfctdef.jump.trlpadding    = automatically determined based on the filter padding (cfg.padding)
%   cfg.artfctdef.jump.artpadding    = automatically determined based on the filter padding (cfg.padding)
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
% cfg.artfctdef.zvalue.sgn in ARTIFACT_JUMP as cfg.artfctdef.jump.sgn default = 'MEG', documented
% cfg.artfctdef.zvalue.cutoff in ARTIFACT_JUMP as cfg.artfctdef.jump.cutoff default = 20, documented
% cfg.artfctdef.zvalue.trlpadding in ARTIFACT_JUMP as cfg.artfctdef.jump.trlpadding automatically determined based on the filter padding (cfg.padding), documented
% cfg.artfctdef.zvalue.fltpadding in ARTIFACT_JUMP as cfg.artfctdef.jump.fltpadding default = 0
% cfg.artfctdef.zvalue.artpadding in ARTIFACT_JUMP as cfg.artfctdef.jump.artpadding automatically determined based on the filter padding (cfg.padding), documented
% cfg.artfctdef.zvalue.feedback

% Copyright (c) 2003-2006, Jan-Mathijs Schoffelen & Robert Oostenveld
%
% $Log: artifact_jump.m,v $
% Revision 1.15  2006/11/29 09:06:36  roboos
% renamed all cfg options with "sgn" into "channel", added backward compatibility when required
% updated documentation, mainly in the artifact detection routines
%
% Revision 1.14  2006/04/25 17:06:28  ingnie
% updated documentation
%
% Revision 1.13  2006/04/20 09:58:33  roboos
% updated documentation
%
% Revision 1.12  2006/02/27 17:01:07  roboos
% added tmpcfg.dataset, added try-end around copying of headerfile, datafile and dataset
%
% Revision 1.11  2006/01/17 14:05:27  roboos
% do preproc absdiff instead of rectify in combination with derivative, absdiff ensures the right order
%
% Revision 1.10  2006/01/12 13:51:38  roboos
% completely new implementation, all based upon the same artifact_zvalue code
% all preprocessing is now done consistently and the various paddings have been better defined
% the functions do not have any explicit support any more for non-continuous data
% the old artifact_xxx functions from JM have been renamed to xxx_old
%

% set default rejection parameters
if ~isfield(cfg,'artfctdef'),                      cfg.artfctdef                 = [];              end
if ~isfield(cfg.artfctdef,'jump'),                 cfg.artfctdef.jump            = [];              end
if ~isfield(cfg.artfctdef.jump,'method'),          cfg.artfctdef.jump.method     = 'zvalue';        end

% for backward compatibility
if isfield(cfg.artfctdef.jump,'sgn')
  cfg.artfctdef.jump.channel = cfg.artfctdef.jump.sgn;
  cfg.artfctdef.jump         = rmfield(cfg.artfctdef.jump, 'sgn');
end

if isfield(cfg.artfctdef.jump, 'artifact') 
  fprintf('jump artifact detection has already been done, retaining artifacts\n');
  artifact = cfg.artfctdef.jump.artifact;
  return
end

if strcmp(cfg.artfctdef.jump.method, 'zvalue')
  % the following fields should be supported for backward compatibility
  dum = 0;
  if isfield(cfg.artfctdef.jump,'pretim'),
    dum = max(dum, cfg.artfctdef.jump.pretim);
    cfg.artfctdef.jump = rmfield(cfg.artfctdef.jump,'pretim');
  end
  if isfield(cfg.artfctdef.jump,'psttim'),
    dum = max(dum, cfg.artfctdef.jump.psttim);
    cfg.artfctdef.jump = rmfield(cfg.artfctdef.jump,'psttim');
  end
  if dum
    cfg.artfctdef.jump.artpadding = max(dum);
  end
  if isfield(cfg.artfctdef.jump,'padding'),
    cfg.artfctdef.jump.trlpadding   = cfg.artfctdef.jump.padding;
    cfg.artfctdef.jump = rmfield(cfg.artfctdef.jump,'padding');
  end
  % settings for preprocessing
  if ~isfield(cfg.artfctdef.jump,'medianfilter'),  cfg.artfctdef.jump.medianfilter  = 'yes';        end
  if ~isfield(cfg.artfctdef.jump,'medianfiltord'), cfg.artfctdef.jump.medianfiltord = 9;            end
  if ~isfield(cfg.artfctdef.jump,'absdiff'),       cfg.artfctdef.jump.absdiff       = 'yes';        end  % compute abs(diff(data)), whereas the order of rectify=yes in combination with derivative=yes would be diff(abs(data)) due to the ordering in preproc
  % settings for the zvalue subfunction
  if ~isfield(cfg.artfctdef.jump,'cutoff'),        cfg.artfctdef.jump.cutoff     = 20;              end
  if ~isfield(cfg.artfctdef.jump,'channel'),       cfg.artfctdef.jump.channel    = 'MEG';           end
  if isfield(cfg, 'padding') && cfg.padding~=0
     if ~isfield(cfg.artfctdef.jump,'trlpadding'), cfg.artfctdef.jump.trlpadding = 0.5*cfg.padding; end
     if ~isfield(cfg.artfctdef.jump,'artpadding'), cfg.artfctdef.jump.artpadding = 0.5*cfg.padding; end
     if ~isfield(cfg.artfctdef.jump,'fltpadding'), cfg.artfctdef.jump.fltpadding = 0;               end
  else
     if ~isfield(cfg.artfctdef.jump,'trlpadding'), cfg.artfctdef.jump.trlpadding = 0;               end
     if ~isfield(cfg.artfctdef.jump,'artpadding'), cfg.artfctdef.jump.artpadding = 0;               end
     if ~isfield(cfg.artfctdef.jump,'fltpadding'), cfg.artfctdef.jump.fltpadding = 0;               end
  end
  % construct a temporary configuration that can be passed onto artifact_zvalue
  tmpcfg                  = [];
  tmpcfg.trl              = cfg.trl;
  try, tmpcfg.dataset     = cfg.dataset; end
  try, tmpcfg.datafile    = cfg.datafile; end
  try, tmpcfg.headerfile  = cfg.headerfile; end
  tmpcfg.artfctdef.zvalue = cfg.artfctdef.jump;
  if isfield(cfg, 'datatype')
    tmpcfg.datatype = cfg.datatype;
  end
  % call the zvalue artifact detection function
  [tmpcfg, artifact] = artifact_zvalue(tmpcfg);
  cfg.artfctdef.jump = tmpcfg.artfctdef.zvalue;
else
  error(sprintf('jump artifact detection only works with cfg.method=''zvalue'''));
end

