function [cfg, artifact] = artifact_eog_old(cfg);

% ARTIFACT_EOG_OLD is deprecated, please use ARTIFACT_EOG instead

% Undocumented local options:
% cfg.artfctdef
% cfg.datafile
% cfg.datatype
% cfg.headerfile
% cfg.trl
% cfg.version

% Copyright (C) 2004, F.C. Donders Centre
%
% $Log: artifact_eog_old.m,v $
% Revision 1.3  2006/04/20 09:56:07  roboos
% removed the (outdated and unclear) documentation, now it only
% contains a comment that the function is deprecated and has been
% replaced by a new one
%

% set default rejection parameters for eog artifacts if necessary.
if ~isfield(cfg,'artfctdef'),            cfg.artfctdef               = [] ;      end;
if ~isfield(cfg.artfctdef,'eog'),        cfg.artfctdef.eog           = [];       end;
if ~isfield(cfg.artfctdef.eog,'sgn'),    cfg.artfctdef.eog.sgn       = {'EOG'};  end;
if ~isfield(cfg.artfctdef.eog,'method'), cfg.artfctdef.eog.method    = 'zvalue'; end;

% for backward-compatibility
if isfield(cfg.artfctdef.eog,'pssbnd'),  
   cfg.artfctdef.eog.bpfreq   = cfg.artfctdef.eog.pssbnd; 
   cfg.artfctdef.eog          = rmfield(cfg.artfctdef.eog,'pssbnd');
end;

% method independent settings
if ~isfield(cfg.artfctdef.eog,'pretim'), cfg.artfctdef.eog.pretim    = 0.1;      end;
if ~isfield(cfg.artfctdef.eog,'psttim'), cfg.artfctdef.eog.psttim    = 0.1;      end;
if ~isfield(cfg.artfctdef.eog,'padding'),cfg.artfctdef.eog.padding   = 0.5;      end;
if ~isfield(cfg.artfctdef.eog,'feedback'),cfg.artfctdef.eog.feedback = 'no';     end;

% method dependent settings
if strcmp(cfg.artfctdef.eog.method,'zvalue'),
    if ~isfield(cfg.artfctdef.eog,'cutoff'),   cfg.artfctdef.eog.cutoff   = 4;      end;
    if ~isfield(cfg.artfctdef.eog,'bpfilter'), cfg.artfctdef.eog.bpfilter = 'yes';  end;
    if ~isfield(cfg.artfctdef.eog,'bpfreq'),   cfg.artfctdef.eog.bpfreq   = [1 15]; end;
    if ~isfield(cfg.artfctdef.eog,'bpfiltord'),cfg.artfctdef.eog.bpfiltord= 4;      end;
    if ~isfield(cfg.artfctdef.eog,'bpfilttype'),cfg.artfctdef.eog.bpfilttype= 'but';  end;
    if ~isfield(cfg.artfctdef.eog,'hilbert'),  cfg.artfctdef.eog.hilbert  = 'yes';  end;
elseif strcmp(cfg.artfctdef.eog.method,'percnt'),
    if ~isfield(cfg.artfctdef.eog,'cutoff'),   cfg.artfctdef.eog.cutoff   = 0.995;  end;
    if ~isfield(cfg.artfctdef.eog,'bpfilter'), cfg.artfctdef.eog.bpfilter = 'yes';  end;
    if ~isfield(cfg.artfctdef.eog,'bpfreq'),   cfg.artfctdef.eog.bpfreq   = [1 15]; end;
    if ~isfield(cfg.artfctdef.eog,'bpfiltord'),cfg.artfctdef.eog.bpfiltord= 4;      end;
    if ~isfield(cfg.artfctdef.eog,'bpfilttype'),cfg.artfctdef.eog.bpfilttype= 'but';  end;
    if ~isfield(cfg.artfctdef.eog,'hilbert'),  cfg.artfctdef.eog.hilbert  = 'yes';  end;
elseif strcmp(cfg.artfctdef.eog.method,'abs'),
    if ~isfield(cfg.artfctdef.eog,'cutoff'),   cfg.artfctdef.eog.cutoff   = 50e-6;  end;
    if ~isfield(cfg.artfctdef.eog,'hpfilter'), cfg.artfctdef.eog.bpfilter = 'yes';  end;
    if ~isfield(cfg.artfctdef.eog,'hpfreq'),   cfg.artfctdef.eog.bpfreq   = 1;      end;
    if ~isfield(cfg.artfctdef.eog,'bpfiltord'),cfg.artfctdef.eog.bpfiltord= 4;      end;
    if ~isfield(cfg.artfctdef.eog,'bpfilttype'),cfg.artfctdef.eog.bpfilttype= 'but';  end;
    if ~isfield(cfg.artfctdef.eog,'hilbert'),  cfg.artfctdef.eog.hilbert  = 'yes';  end;
end    

if ~isfield(cfg, 'datatype') || ~strcmp(cfg.datatype, 'continuous')
  % datatype is unknown or not continuous, perform epoch boundary check
  iscontinuous = 0;
else
  % do not perform epoch boundary check, usefull for pseudo-continuous data
  iscontinuous = strcmp(cfg.datatype, 'continuous');
end

artfctdef    = cfg.artfctdef.eog;
cfg          = dataset2files(cfg);
hdr          = read_fcdc_header(cfg.headerfile); 
padsmp       = round(artfctdef.padding*hdr.Fs);
trl          = cfg.trl;
trl(:,1)     = trl(:,1) - padsmp; % pad the trial with padsmp-samples, in order to detect
trl(:,2)     = trl(:,2) + padsmp; % artifacts at the edges of the relevant trials.
numtrl       = size(trl,1);
trllength    = trl(:,2) - trl(:,1) + 1;
artfctdef.trl= trl;
artfctdef.sgn= channelselection(artfctdef.sgn,hdr.label);
sgnind       = match_str(hdr.label, artfctdef.sgn);
numeogsgn    = length(sgnind);

if numeogsgn<1
  error('no EOG channels selected');
end

if ~isfield(artfctdef, 'fltpadding')
  if iscontinuous
    if ~isfield(artfctdef, 'bpfilttype') || ~isfield(artfctdef, 'hpfilttype'),
       fltpadding = 0;
    elseif strcmp(artfctdef.bpfilttype, 'but') || strcmp(artfctdef.hpfilttype, 'but'),
      fltpadding = 20*artfctdef.bpfiltord;
    elseif strcmp(artfctdef.bpfilttype, 'fir') || strcmp(artfctdef.hpfilttype, 'fir'),
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

% initialize some variables
fltdatavg    = zeros(numeogsgn,1);
fltdatvar    = zeros(numeogsgn,1);
fltdatstd    = zeros(numeogsgn,1);
cumpernumsmp = 0;
cutfltdatstr = cell(numtrl,1);
artfctchn    = zeros(numtrl,1);
rejectall    = zeros(1,trl(end,2));
z_tdata      = zeros(1,trl(end,2));

for trllop = 1:numtrl  
     fprintf('scanning for eog artifacts in trial %d from %d\n', trllop, numtrl);
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

for trllop = 1:numtrl
    if strcmp(artfctdef.method,'zvalue'),
        dummy = ((cutfltdatstr{trllop} - repmat(fltdatavg,1,trllength(trllop))) ./ ...
                repmat(fltdatstd,1,trllength(trllop)));
        [dummy, index] = sort(dummy,1);
        [dum indx] = max(dummy(end,:));
        artfctchn(trllop)=index(end,indx);
        z_tdata(1,trl(trllop,1):trl(trllop,2)) = sum(dummy,1)./sqrt(size(dummy,1));
    elseif strcmp(artfctdef.method,'abs'),
        dummy = cutfltdatstr{trllop};
        [dummy,index] = sort(dummy);
        [dum indx] = max(dummy(end,:));
        artfctchn(trllop)=index(end,indx);
        z_tdata(1,trl(trllop,1):trl(trllop,2)) = sum(dummy,1)./sqrt(size(dummy,1));
    end;    
end;

if strcmp(artfctdef.method,'percnt'),
    dummy = zeros(numeogsgn,trl(end,2));
    for trllop = 1:numtrl
        dummy(:,trl(trllop,1):trl(trllop,2)) = cutfltdatstr{trllop};
        clear cutfltdatstr{trllop};    
    end;
    
    for sgnlop = 1:numeogsgn
            [sd,id]         = sort(dummy(sgnlop,:),2);
            sel             = find(sd);
            dummy(sgnlop,id(sel)) = [1:length(sel)]./length(sel);                
    end
    [dummy, index] = sort(dummy,1);
    
    for trllop = 1:numtrl
        dum = dummy(:,trl(trllop,1):trl(trllop,2));
        [dum,index] = sort(dum,1);
        artfctchn(trllop)=index(end,find(dum(end,:)==max(dum(end,:))));
        z_tdata(1,trl(trllop,1):trl(trllop,2)) = dum(end,:);
    end
end

clear cutfltdatstr;

% update the configuration prior to calling artifact_feedback
cfg.artfctdef.eog = artfctdef;

if (strcmp(cfg.artfctdef.eog.feedback,'yes'))
    [artfctdef,artifact] = artifact_feedback_old(cfg, z_tdata, artfctchn,'EOGraw');
elseif (strcmp(cfg.artfctdef.eog.feedback,'flt'))
    [artfctdef,artifact] = artifact_feedback_old(cfg, z_tdata, artfctchn,'EOGflt');
else    
    [artfctdef,artifact] = artifact_feedback_old(cfg, z_tdata, artfctchn,'EOG');   
end;

% remember the details that were used here
cfg.artfctdef.eog          = artfctdef;
cfg.artfctdef.eog.artifact = artifact;

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: artifact_eog_old.m,v 1.3 2006/04/20 09:56:07 roboos Exp $';

