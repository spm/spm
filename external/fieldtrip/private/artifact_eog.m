function [cfg, artifact] = artifact_eog(cfg,data)

% ARTIFACT_EOG reads the data segments of interest from file and
% identifies EOG artifacts.
%
% Use as
%   [cfg, artifact] = artifact_eog(cfg)
% or
%   [cfg, artifact] = artifact_eog(cfg, data)
%
% If you are calling ARTIFACT_EOG with only the configuration as first
% input argument and the data still has to be read from file, you should
% specify:
%   cfg.dataset    = string, name of the dataset
% or, instead of cfg.dataset
%   cfg.datafile   = string
%   cfg.headerfile = string
%
% If you are calling ARTIFACT_EOG with also the second input argument
% "data", then that should contain data that was already read from file in
% a call to PREPROCESSING.
%
% Always specify:
%   cfg.trl        = structure that defines the data segments of interest. See DEFINETRIAL
%
% The data is preprocessed (again) with the following configuration parameters,
% which are optimal for identifying EOG artifacts:
%   cfg.artfctdef.eog.bpfilter   = 'yes'
%   cfg.artfctdef.eog.bpfilttype = 'but'
%   cfg.artfctdef.eog.bpfreq     = [1 15]
%   cfg.artfctdef.eog.bpfiltord  = 4
%   cfg.artfctdef.eog.hilbert    = 'yes'
%
% Artifacts are identified by means of thresholding the z-transformed value
% of the preprocessed data.
%   cfg.artfctdef.eog.channel      = Nx1 cell-array with selection of channels, see CHANNELSELECTION for details
%   cfg.artfctdef.eog.cutoff       = 4       z-value at which to threshold
%   cfg.artfctdef.eog.trlpadding   = 0.5
%   cfg.artfctdef.eog.fltpadding   = 0.1
%   cfg.artfctdef.eog.artpadding   = 0.1
%
% The output argument "artifact" is a Nx2 matrix comparable to the
% "trl" matrix of DEFINETRIAL. The first column of which specifying the
% beginsamples of an artifact period, the second column contains the
% endsamples of the artifactperiods.
%
% See also ARTIFACT_ZVALUE, REJECTARTIFACT

% Undocumented local options
% cfg.datatype
% cfg.method

% Copyright (c) 2003-2006, Jan-Mathijs Schoffelen & Robert Oostenveld
%
% $Log: artifact_eog.m,v $
% Revision 1.29  2008/10/10 15:00:10  estmee
% Repaired determining artifacts with cfg as only input argument.
%
% Revision 1.28  2008/10/07 16:11:08  estmee
% Added data as second input argument to artifact_eog itself and to the way it calls artifact_zvalue.
%
% Revision 1.27  2008/10/07 08:58:51  roboos
% committed the changes that Esther made recently, related to the support of data as input argument to the artifact detection functions. I hope that this does not break the functions too seriously.
%
% Revision 1.26  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.25  2006/11/29 09:06:36  roboos
% renamed all cfg options with "sgn" into "channel", added backward compatibility when required
% updated documentation, mainly in the artifact detection routines
%
% Revision 1.24  2006/04/25 17:06:28  ingnie
% updated documentation
%
% Revision 1.23  2006/04/20 09:58:33  roboos
% updated documentation
%
% Revision 1.22  2006/02/27 17:01:07  roboos
% added tmpcfg.dataset, added try-end around copying of headerfile, datafile and dataset
%
% Revision 1.21  2006/01/12 13:51:38  roboos
% completely new implementation, all based upon the same artifact_zvalue code
% all preprocessing is now done consistently and the various paddings have been better defined
% the functions do not have any explicit support any more for non-continuous data
% the old artifact_xxx functions from JM have been renamed to xxx_old
%

fieldtripdefs

% set default rejection parameters
if ~isfield(cfg,'artfctdef'),                  cfg.artfctdef                 = [];       end
if ~isfield(cfg.artfctdef,'eog'),              cfg.artfctdef.eog             = [];       end
if ~isfield(cfg.artfctdef.eog,'method'),       cfg.artfctdef.eog.method      = 'zvalue'; end

% for backward compatibility
if isfield(cfg.artfctdef.eog,'sgn')
  cfg.artfctdef.eog.channel = cfg.artfctdef.eog.sgn;
  cfg.artfctdef.eog         = rmfield(cfg.artfctdef.eog, 'sgn');
end

if isfield(cfg.artfctdef.eog, 'artifact')
  fprintf('eog artifact detection has already been done, retaining artifacts\n');
  artifact = cfg.artfctdef.eog.artifact;
  return
end

if strcmp(cfg.artfctdef.eog.method, 'zvalue')
  % the following fields should be supported for backward compatibility
  if isfield(cfg.artfctdef.eog,'pssbnd'),
    cfg.artfctdef.eog.bpfreq   = cfg.artfctdef.eog.pssbnd;
    cfg.artfctdef.eog.bpfilter = 'yes';
    cfg.artfctdef.eog = rmfield(cfg.artfctdef.eog,'pssbnd');
  end;
  dum = 0;
  if isfield(cfg.artfctdef.eog,'pretim'),
    dum = max(dum, cfg.artfctdef.eog.pretim);
    cfg.artfctdef.eog = rmfield(cfg.artfctdef.eog,'pretim');
  end
  if isfield(cfg.artfctdef.eog,'psttim'),
    dum = max(dum, cfg.artfctdef.eog.psttim);
    cfg.artfctdef.eog = rmfield(cfg.artfctdef.eog,'psttim');
  end
  if dum
    cfg.artfctdef.eog.artpadding = max(dum);
  end
  if isfield(cfg.artfctdef.eog,'padding'),
    cfg.artfctdef.eog.trlpadding   = cfg.artfctdef.eog.padding;
    cfg.artfctdef.eog = rmfield(cfg.artfctdef.eog,'padding');
  end
  % settings for preprocessing
  if ~isfield(cfg.artfctdef.eog,'bpfilter'),   cfg.artfctdef.eog.bpfilter   = 'yes';     end
  if ~isfield(cfg.artfctdef.eog,'bpfilttype'), cfg.artfctdef.eog.bpfilttype = 'but';     end
  if ~isfield(cfg.artfctdef.eog,'bpfreq'),     cfg.artfctdef.eog.bpfreq     = [1 15];    end
  if ~isfield(cfg.artfctdef.eog,'bpfiltord'),  cfg.artfctdef.eog.bpfiltord  = 4;         end
  if ~isfield(cfg.artfctdef.eog,'hilbert'),    cfg.artfctdef.eog.hilbert    = 'yes';     end
  % settings for the zvalue subfunction
  if ~isfield(cfg.artfctdef.eog,'channel'),    cfg.artfctdef.eog.channel     = 'EOG';    end
  if ~isfield(cfg.artfctdef.eog,'trlpadding'), cfg.artfctdef.eog.trlpadding  = 0.5;      end
  if ~isfield(cfg.artfctdef.eog,'artpadding'), cfg.artfctdef.eog.artpadding  = 0.1;      end
  if ~isfield(cfg.artfctdef.eog,'fltpadding'), cfg.artfctdef.eog.fltpadding  = 0.1;      end
  if ~isfield(cfg.artfctdef.eog,'cutoff'),     cfg.artfctdef.eog.cutoff      = 4;        end
  % construct a temporary configuration that can be passed onto artifact_zvalue
  tmpcfg                  = [];
  tmpcfg.trl              = cfg.trl;
  try, tmpcfg.dataset     = cfg.dataset; end
  try, tmpcfg.datafile    = cfg.datafile; end
  try, tmpcfg.headerfile  = cfg.headerfile; end
  tmpcfg.artfctdef.zvalue = cfg.artfctdef.eog;
  if isfield(cfg, 'datatype')
    tmpcfg.datatype = cfg.datatype;
  end
  % call the zvalue artifact detection function
  if nargin ==1
    [tmpcfg, artifact] = artifact_zvalue(tmpcfg);  
  elseif nargin ==2
    [tmpcfg, artifact] = artifact_zvalue(tmpcfg,data);      
  end
cfg.artfctdef.eog  = tmpcfg.artfctdef.zvalue;
else
  error(sprintf('EOG artifact detection only works with cfg.method=''zvalue'''));
end

