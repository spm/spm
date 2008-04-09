function [cfg, artifact] = artifact_muscle(cfg);

% ARTIFACT_MUSCLE_OLD is deprecated, please use ARTIFACT_MUSCLE instead

% Undocumented local options:
% cfg.artfctdef
% cfg.datafile
% cfg.datatype
% cfg.headerfile
% cfg.trl
% cfg.version
% Copyright (c) 2003, F.C. Donders Centre
%
% $Log: artifact_muscle_old.m,v $
% Revision 1.3  2006/04/20 09:56:07  roboos
% removed the (outdated and unclear) documentation, now it only
% contains a comment that the function is deprecated and has been
% replaced by a new one
%

%set default rejection parameters for muscle artifacts if necessary.          
if ~isfield(cfg,'artfctdef'),               cfg.artfctdef                  = [];      end;
if ~isfield(cfg.artfctdef,'muscle'),        cfg.artfctdef.muscle           = [];      end;
if ~isfield(cfg.artfctdef.muscle,'sgn'),    cfg.artfctdef.muscle.sgn       = 'MEG';   end;
if ~isfield(cfg.artfctdef.muscle,'pretim'), cfg.artfctdef.muscle.pretim    = 0.05;    end;
if ~isfield(cfg.artfctdef.muscle,'psttim'), cfg.artfctdef.muscle.psttim    = 0.05;    end;
if ~isfield(cfg.artfctdef.muscle,'padding'),cfg.artfctdef.muscle.padding   = 0.1;     end;
if ~isfield(cfg.artfctdef.muscle,'method'), cfg.artfctdef.muscle.method    = 'zvalue';end;
if ~isfield(cfg.artfctdef.muscle,'feedback'),cfg.artfctdef.muscle.feedback = 'no';    end;

% for backward-compatibility
if isfield(cfg.artfctdef.muscle,'pssbnd'),  
   cfg.artfctdef.muscle.bpfreq   = cfg.artfctdef.muscle.pssbnd; 
   cfg.artfctdef.muscle          = rmfield(cfg.artfctdef.muscle,'pssbnd');
end;

% method dependent settings
if strcmp(cfg.artfctdef.muscle.method,'zvalue'),
    if ~isfield(cfg.artfctdef.muscle,'bpfilter'),cfg.artfctdef.muscle.bpfilter = 'yes';   end;
    if ~isfield(cfg.artfctdef.muscle,'bpfreq'), cfg.artfctdef.muscle.bpfreq = [110 140];  end;
    if ~isfield(cfg.artfctdef.muscle,'bpfiltord'), cfg.artfctdef.muscle.bpfiltord = 10;         end;
    if ~isfield(cfg.artfctdef.muscle,'bpfilttype'),cfg.artfctdef.muscle.bpfilttype= 'but';      end;
    if ~isfield(cfg.artfctdef.muscle,'cutoff'), cfg.artfctdef.muscle.cutoff = 4;          end;
    if ~isfield(cfg.artfctdef.muscle,'hilbert'),cfg.artfctdef.muscle.hilbert = 'yes';     end;      
else 
    error('muscle-artifact rejection only works with method: zvalue\n');
end

if ~isfield(cfg, 'datatype') || ~strcmp(cfg.datatype, 'continuous')
  % datatype is unknown or not continuous, perform epoch boundary check
  iscontinuous = 0;
else
  % do not perform epoch boundary check, usefull for pseudo-continuous data
  iscontinuous = strcmp(cfg.datatype, 'continuous');
end

artfctdef    = cfg.artfctdef.muscle;
cfg          = dataset2files(cfg);
hdr          = read_fcdc_header(cfg.headerfile);
padsmp       = round(artfctdef.padding*hdr.Fs);
trl          = cfg.trl;
trl(:,1)     = trl(:,1) - padsmp; %pad the trial with padsmp-samples, in order to detect
trl(:,2)     = trl(:,2) + padsmp; %artifacts at the edges of the relevant trials.
numtrl       = size(trl,1);
trllength    = trl(:,2) - trl(:,1) + 1;%compute number of samples per trial
artfctdef.trl= trl;
artfctdef.sgn= channelselection(artfctdef.sgn,hdr.label);
sgnind       = match_str(hdr.label, artfctdef.sgn);
nummusclesgn = length(sgnind);

if ~isfield(artfctdef, 'fltpadding')
  if iscontinuous
    if ~isfield(artfctdef, 'bpfilttype')
      fltpadding = 0;
    elseif strcmp(artfctdef.bpfilttype, 'but')
      fltpadding = 20*artfctdef.bpfiltord;
    elseif strcmp(artfctdef.bpfilttype, 'fir')
      fltpadding = artfctdef.bpfiltord;
    else
      warning('unknown filter type, cannot determine filter padding');
      fltpadding = 0;
    end
  else
    warning('default is not to apply filter padding on trial based data');
    fltpadding = 0;
  end
else
  fltpadding = round(artfctdef.fltpadding*hdr.Fs); 
end

fltdatavg    = zeros(nummusclesgn,1);
fltdatvar    = zeros(nummusclesgn,1);
fltdatstd    = zeros(nummusclesgn,1);
cumpernumsmp = 0;
cutfltdatstr = cell(numtrl,1);
cutdatstr = cell(numtrl,1);
artfctchn    = zeros(numtrl,1);
rejectall    = zeros(1,trl(end,2));
z_tdata      = zeros(1,trl(end,2));

%loop over all trials and calculate mean and std
for trllop = 1:numtrl  
     fprintf('scanning for muscle artifacts in trial %d from %d\n', trllop, numtrl);
     dat = read_fcdc_data(cfg.datafile, hdr, trl(trllop,1)-fltpadding, ...
         trl(trllop,2)+fltpadding, sgnind, iscontinuous);%read data
     fltdat = preproc(dat, artfctdef.sgn, hdr.Fs, artfctdef, [], fltpadding, fltpadding);
     cutfltdatstr{trllop} = fltdat; 
     fltdatavg = fltdatavg + mean(fltdat,2) .* trllength(trllop); 
     fltdatvar = fltdatvar + (std(fltdat,[],2)).^2 .* trllength(trllop); 
     cumpernumsmp = cumpernumsmp + trllength(trllop);
end

fltdatavg = fltdatavg ./ cumpernumsmp;
fltdatstd = sqrt(fltdatvar ./ cumpernumsmp);

%compute z-values for each channel and calculate cumulative z-scores for each sample
%keep track of the channel with the highest individual z-value in each trial for displaying purposes.
for trllop = 1:numtrl
    dummy = ((cutfltdatstr{trllop} - repmat(fltdatavg,1,trllength(trllop))) ./ ...
        repmat(fltdatstd,1,trllength(trllop)));
    [dummy,index]=sort(dummy,1);
    [dum indx] = max(dummy(end,:));
    artfctchn(trllop)=index(end,indx);    
    z_tdata(1,trl(trllop,1):trl(trllop,2)) = sum(dummy,1)./sqrt(size(dummy,1));
end;

clear cutfltdatstr;

% update the configuration prior to calling artifact_feedback
cfg.artfctdef.muscle = artfctdef;

[artfctdef, artifact] = artifact_feedback_old(cfg, z_tdata, artfctchn, 'muscle');

% remember the details that were used here
cfg.artfctdef.muscle          = artfctdef;
cfg.artfctdef.muscle.artifact = artifact;

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: artifact_muscle_old.m,v 1.3 2006/04/20 09:56:07 roboos Exp $';

