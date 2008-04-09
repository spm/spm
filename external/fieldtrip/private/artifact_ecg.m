function [cfg, artifact] = artifact_ecg(cfg)

% ARTIFACT_ECG performs a peak-detection on the ECG-channel. The
% heart activity can be seen in the MEG data as an MCG artifact and
% can be removed using independent component analysis.
% 
% Use as
%   [cfg, artifact] = artifact_ecg(cfg)
% where the configuration should contain
%   cfg.artfctdef.ecg.channel = Nx1 cell-array with selection of channels, see CHANNELSELECTION for details
%   cfg.artfctdef.ecg.pretim  = 0.05; pre-artifact rejection-interval in seconds
%   cfg.artfctdef.ecg.psttim  = 0.3;  post-artifact rejection-interval in seconds
%   cfg.artfctdef.ecg.method  = 'zvalue'; peak-detection method
%   cfg.artfctdef.ecg.cutoff  = 3; peak-threshold
%
% The output artifact variable is an Nx2-matrix, containing the 
% begin and end samples of the QRST-complexes in the ECG.
%
% See also REJECTARTIFACT

% Undocumented local options:
% cfg.datatype
% cfg.trl
% cfg.version
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

% Copyright (c) 2005, Jan-Mathijs Schoffelen
%
% $Log: artifact_ecg.m,v $
% Revision 1.12  2007/07/31 08:40:59  jansch
% added option cfg.arfctdef.ecg.mindist to specify the minimal distance between
% the peaks
%
% Revision 1.11  2006/11/29 09:06:36  roboos
% renamed all cfg options with "sgn" into "channel", added backward compatibility when required
% updated documentation, mainly in the artifact detection routines
%
% Revision 1.10  2006/06/14 12:43:49  roboos
% removed the documentation for cfg.lnfilttype, since that option is not supported by preproc
%
% Revision 1.9  2006/04/20 09:58:33  roboos
% updated documentation
%
% Revision 1.8  2006/04/12 08:38:01  ingnie
% updated documentation
%
% Revision 1.7  2006/01/31 13:49:29  jansch
% included dataset2files to ensure the presence of cfg.headerfile or cfg.datafile
% whenever needed
%
% Revision 1.6  2006/01/30 14:07:01  jansch
% added new ntrl in line 135 since the number of trials could have changed
%
% Revision 1.5  2005/12/20 08:36:47  roboos
% add the artifact Nx2 matrix to the output configuration
% changed some indentation and white space, renamed a few variables
%
% Revision 1.4  2005/09/29 14:56:32  jansch
% first implementation
%

% set default rejection parameters for eog artifacts if necessary.
if ~isfield(cfg,'artfctdef'),            cfg.artfctdef               = [];            end
if ~isfield(cfg.artfctdef,'ecg'),        cfg.artfctdef.ecg           = [];            end
if ~isfield(cfg.artfctdef.ecg,'channel'),cfg.artfctdef.ecg.channel   = {'ECG'};       end
if ~isfield(cfg.artfctdef.ecg,'method'), cfg.artfctdef.ecg.method    = 'zvalue';      end
if ~isfield(cfg.artfctdef.ecg,'cutoff'), cfg.artfctdef.ecg.cutoff    = 3;             end
if ~isfield(cfg.artfctdef.ecg,'padding'),cfg.artfctdef.ecg.padding   = 0.5;           end
if ~isfield(cfg.artfctdef.ecg,'inspect'),cfg.artfctdef.ecg.inspect   = {'MLT' 'MRT'}; end
if ~isfield(cfg.artfctdef.ecg,'pretim'), cfg.artfctdef.ecg.pretim    = 0.05;          end
if ~isfield(cfg.artfctdef.ecg,'psttim'), cfg.artfctdef.ecg.psttim    = 0.3;           end
if ~isfield(cfg.artfctdef.ecg,'mindist'), cfg.artfctdef.ecg.mindist  = 0.5;           end

% for backward compatibility
if isfield(cfg.artfctdef.ecg,'sgn')
  cfg.artfctdef.ecg.channel = cfg.artfctdef.ecg.sgn;
  cfg.artfctdef.ecg         = rmfield(cfg.artfctdef.ecg, 'sgn');
end

if ~strcmp(cfg.artfctdef.ecg.method, 'zvalue'),
  error('this method is not applicable');
end

if ~isfield(cfg, 'datatype') || ~strcmp(cfg.datatype, 'continuous')
  % datatype is unknown or not continuous, perform epoch boundary check
  iscontinuous = 0;
  error('not implemented yet for discontinuous data');
else
  % do not perform epoch boundary check, usefull for pseudo-continuous data
  iscontinuous = strcmp(cfg.datatype, 'continuous');
end

artfctdef     = cfg.artfctdef.ecg;
cfg           = dataset2files(cfg);
hdr           = read_fcdc_header(cfg.headerfile);
padsmp        = round(artfctdef.padding*hdr.Fs);
trl           = cfg.trl;
ntrl          = size(trl,1);
artfctdef.trl = trl;
artfctdef.channel = channelselection(artfctdef.channel, hdr.label);
artfctdef.blc = 'yes';
sgnind        = match_str(hdr.label, artfctdef.channel);
numecgsgn     = length(sgnind);
fltpadding    = 0;

if numecgsgn<1
  error('no ECG channels selected');
elseif numecgsgn>1
  error('only one ECG channel can be selected');
end

% read in the ecg-channel and do blc and squaring
for j = 1:ntrl
  ecg{j} = read_fcdc_data(cfg.datafile, hdr, trl(j,1), trl(j,2), sgnind, iscontinuous);
  ecg{j} = preproc(ecg{j}, artfctdef.channel, hdr.Fs, artfctdef, [], fltpadding, fltpadding);
  ecg{j} = ecg{j}.^2;
end 

tmp   = cell2mat(ecg);
stmp  =  std(tmp, 0, 2);
mtmp  = mean(tmp, 2);
Nsmp  = max(trl(:,2));
trace = zeros(1,Nsmp);

% standardise the ecg
for j = 1:ntrl
  trace(trl(j,1):trl(j,2)) = (ecg{j}-mtmp)./stmp;
end

accept = 0;
while accept == 0,
  h = figure;  
  plot(trace);
  hold on; 
  plot([1 Nsmp], [artfctdef.cutoff artfctdef.cutoff], 'r:');
  hold off;
  xlabel('samples');
  ylabel('zscore');

  fprintf(['\ncurrent  ',artfctdef.method,' threshold = %1.3f'], artfctdef.cutoff);
  response = input('\nkeep the current value (y/n) ?\n','s');
  switch response
  case 'n'
      oldcutoff = artfctdef.cutoff;
      artfctdef.cutoff = input('\nenter new value \n');
  case 'y'
      oldcutoff = artfctdef.cutoff;
      accept = 1;
  otherwise
      warning('unrecognised response, assuming no');
      oldcutoff = artfctdef.cutoff;
      artfctdef.cutoff = input('\nenter new value \n');
  end;
  close
end

% detect peaks which are at least half a second apart and store
% the indices of the qrs-complexes in the artifact-configuration
mindist       = round(cfg.artfctdef.ecg.mindist.*hdr.Fs);
[pindx, pval] = peakdetect2(trace, artfctdef.cutoff, mindist);
%sel           = find(standardise(pval,2)<2);
%pindx         = pindx(sel);
%pval          = pval(sel);
artfctdef.qrs = pindx;

%---------------------------------------
% create trials for qrs-triggered average
trl = [];
trl(:,1) = pindx(:) - artfctdef.padding*(hdr.Fs)  ;
trl(:,2) = pindx(:) + artfctdef.padding*(hdr.Fs)-1;
trl(:,3) = -artfctdef.padding*(hdr.Fs);
%------------

% ---------------------
% qrs-triggered average
% FIXME, at present this only works for continuous data: the assumption can
% be made that all trials are equally long.
sgn    = channelselection(artfctdef.inspect, hdr.label);
megind = match_str(hdr.label, sgn);
sgnind = [megind(:); sgnind];
dat    = zeros(length(sgnind), trl(1,2)-trl(1,1)+1);
ntrl   = size(trl,1);

if ~isempty(sgnind)
  for j = 1:ntrl
    fprintf('reading and preprocessing trial %d of %d\n', j, ntrl);
    dum = read_fcdc_data(cfg.datafile, hdr, trl(j,1), trl(j,2), sgnind, iscontinuous);
    dat = dat + blc(dum);
  end 
end

dat  = dat./ntrl;
time = offset2time(trl(1,3), hdr.Fs, size(dat,2));
tmp  = dat(1:end-1,:);
mdat = max(abs(tmp(:)));

acceptpre = 0;
acceptpst = 0;
while acceptpre == 0 || acceptpst == 0,
  h = figure;
  subplot(2,1,1); plot(time, dat(end, :));
  subplot(2,1,2);
  axis([time(1) time(end) -1.1*mdat 1.1*mdat]);
  xpos   = -artfctdef.pretim;
  ypos   = -1.05*mdat;
  width  = artfctdef.pretim + artfctdef.psttim;
  height = 2.1*mdat;
  rectangle('Position', [xpos ypos width height], 'FaceColor', 'r');
  hold on; plot(time, dat(1:end-1, :), 'b');
  
  if acceptpre == 0, 
    fprintf(['\ncurrent pre-peak interval = %1.3f'], artfctdef.pretim);
    response = input('\nkeep the current value (y/n) ?\n','s');
    switch response
    case 'n'
        oldpretim = artfctdef.pretim;
        artfctdef.pretim = input('\nenter new value \n');
    case 'y'
        oldpretim = artfctdef.pretim;
        acceptpre = 1;
    otherwise
        warning('unrecognised response, assuming no');
        oldpretim = artfctdef.pretim;
    end     
  end
  if acceptpst == 0 && acceptpre == 1,
    fprintf(['\ncurrent post-peak interval = %1.3f'], artfctdef.psttim);
    response = input('\nkeep the current value (y/n) ?\n','s');
    switch response
    case 'n'
        oldpsttim = artfctdef.psttim;
        artfctdef.psttim = input('\nenter new value \n');
    case 'y'
        oldpsttim = artfctdef.psttim;
        acceptpst = 1;
    otherwise
        warning('unrecognised response, assuming no');
        oldpsttim = artfctdef.psttim;
    end
  end
  close
end

artifact(:,1) = trl(:,1) - trl(:,3) - round(artfctdef.pretim*hdr.Fs);
artifact(:,2) = trl(:,1) - trl(:,3) + round(artfctdef.psttim*hdr.Fs);

% remember the details that were used here
cfg.artfctdef.ecg          = artfctdef;
cfg.artfctdef.ecg.artifact = artifact;

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: artifact_ecg.m,v 1.12 2007/07/31 08:40:59 jansch Exp $';

