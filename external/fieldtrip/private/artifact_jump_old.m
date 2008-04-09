function [cfg, artifact] = artifact_jump_old(cfg);

% ARTIFACT_JUMP_OLD is deprecated, please use ARTIFACT_JUMP instead

% Undocumented local options:
% cfg.artfctdef
% cfg.datafile
% cfg.datatype
% cfg.headerfile
% cfg.padding
% cfg.trl
% cfg.version

% Copyright (C) 2004, F.C. Donders Centre
%
% $Log: artifact_jump_old.m,v $
% Revision 1.3  2006/04/20 09:56:07  roboos
% removed the (outdated and unclear) documentation, now it only
% contains a comment that the function is deprecated and has been
% replaced by a new one
%

% set default rejection parameters for jump artifacts if necessary.          
if ~isfield(cfg,'artfctdef'),               cfg.artfctdef               = [];              end;
if ~isfield(cfg.artfctdef,'jump'),          cfg.artfctdef.jump          = [];              end;
if ~isfield(cfg.artfctdef.jump,'sgn'),      cfg.artfctdef.jump.sgn      = 'MEG';           end;
if ~isfield(cfg.artfctdef.jump,'method'),   cfg.artfctdef.jump.method   = 'zvalue';        end;
if ~isfield(cfg.artfctdef.jump,'feedback'), cfg.artfctdef.jump.feedback = 'no';            end;
if ~isfield(cfg.artfctdef.jump,'medianfilter'), cfg.artfctdef.jump.medianfilter  = 'yes';  end;
if ~isfield(cfg.artfctdef.jump,'medianfiltord'),cfg.artfctdef.jump.medianfiltord = 9;      end;

% padding dependent settings
if ~isfield(cfg,'padding'),	             cfg.padding = 0;                              end;
if cfg.padding~=0 
   if ~isfield(cfg.artfctdef.jump,'pretim'), cfg.artfctdef.jump.pretim  = cfg.padding;     end;
   if ~isfield(cfg.artfctdef.jump,'psttim'), cfg.artfctdef.jump.psttim  = cfg.padding;     end;
   if ~isfield(cfg.artfctdef.jump,'padding'),cfg.artfctdef.jump.padding = 0.5*cfg.padding; end;
else
   if ~isfield(cfg.artfctdef.jump,'pretim'), cfg.artfctdef.jump.pretim  = 0;               end;
   if ~isfield(cfg.artfctdef.jump,'psttim'), cfg.artfctdef.jump.psttim  = 0;               end;
   if ~isfield(cfg.artfctdef.jump,'padding'),cfg.artfctdef.jump.padding = 0;               end;
end;

% method dependent settings
artfctdef.fltord = 0;
if strcmp(cfg.artfctdef.jump.method,'zvalue'),
    if ~isfield(cfg.artfctdef.jump,'cutoff'), cfg.artfctdef.jump.cutoff    = 20;           end;
elseif strcmp(cfg.artfctdef.jump.method,'abs'),
    if ~isfield(cfg.artfctdef.jump,'cutoff'), cfg.artfctdef.jump.cutoff    = 5e-11;        end;
    if ~isfield(cfg.artfctdef.jump,'blc'),    cfg.artfctdef.jump.blc       = 'yes';        end;
end

if ~isfield(cfg, 'datatype') || ~strcmp(cfg.datatype, 'continuous')
  % datatype is unknown or not continuous, perform epoch boundary check
  iscontinuous = 0;
else
  % do not perform epoch boundary check, usefull for pseudo-continuous data
  iscontinuous = strcmp(cfg.datatype, 'continuous');
end

fprintf('\nscanning for jump artifacts');

artfctdef     = cfg.artfctdef.jump;
cfg           = dataset2files(cfg);
hdr           = read_fcdc_header(cfg.headerfile);
trl           = cfg.trl;
numtrl        = size(trl,1);
padding       = round(artfctdef.padding*hdr.Fs);
trl(:,1)      = trl(:,1)-padding;
trl(:,2)      = trl(:,2)+padding;
if ~isempty(find(trl(:,1)<1)),
  fprintf('\nWARNING: required padding goes beyond begin of dataset');   
  trl(find(trl(:,1)<1),1) = 1;
end  
nsmp          = trl(:,2) - trl(:,1) + 1;%compute number of samples per trial
artfctdef.sgn = channelselection(artfctdef.sgn,hdr.label);
sgnind        = match_str(hdr.label, artfctdef.sgn);
artfctdef.trl = trl;
cutfltdatstr  = cell(numtrl,1);
artfctchn     = zeros(numtrl,1);
z_tdata       = zeros(1,trl(end,2));

if ~isfield(artfctdef, 'fltpadding')
  fltpadding = 0;
else
  fltpadding = round(artfctdef.fltpadding*hdr.Fs); 
end

% loop over all signals on which jump-artifact detection is required.
for sgnlop = 1:length(sgnind)
   fprintf('\nreading in signal %d from %d', sgnlop, length(sgnind)); 
   fltdatavg    = 0;
   fltdatstd    = 0;
   fltdatvar    = 0;
   cumpernumsmp = 0;
   % loop over the trials and calculate mean and std
   for trllop = 1:numtrl 
       dat = read_fcdc_data(cfg.datafile, hdr, trl(trllop,1)-fltpadding, ...
            trl(trllop,2)+fltpadding, sgnind(sgnlop), iscontinuous);%read data
       if strcmp(artfctdef.method,'zvalue'),
           fltdat = preproc(dat, artfctdef.sgn, hdr.Fs, artfctdef, [], fltpadding, fltpadding);
           fltdat = [abs(diff(fltdat,1,2)) 0];
       elseif strcmp(artfctdef.method,'abs'),
            fltdat = preproc(dat, artfctdef.sgn, hdr.Fs, artfctdef, [], fltpadding, fltpadding);
           fltdat = fltdat - (min(fltdat));
           fltdat = fltdat(1:end-1);
       end;
       cutfltdatstr{trllop} = fltdat; 
       fltdatavg = fltdatavg + mean(fltdat,2) .* nsmp(trllop); 
       fltdatvar = fltdatvar + (std(fltdat,[],2)).^2 .* nsmp(trllop); 
       cumpernumsmp = cumpernumsmp + nsmp(trllop);
   end

    fltdatavg = fltdatavg ./ cumpernumsmp;
    fltdatstd = sqrt(fltdatvar ./ cumpernumsmp);
	
    % select for each sample the highest z-value, and keep track of the highest z-value for the
    % respective channel in each trial
    for trllop = 1:numtrl
        if strcmp(artfctdef.method,'zvalue'), 
            z_tdat = (cutfltdatstr{trllop} - repmat(fltdatavg,[1 nsmp(trllop)])) ./ ...
                repmat(fltdatstd,[1 nsmp(trllop)]);
        elseif strcmp(artfctdef.method,'abs'),
            z_tdat = cutfltdatstr{trllop};    
        end;
        dummy = sort([z_tdat;z_tdata(1,trl(trllop,1):trl(trllop,2))]);
        z_tdata(1,trl(trllop,1):trl(trllop,2)) = dummy(end,:); 
        maxval(trllop,sgnlop) = max(z_tdat);
    end;
end;

% compute for each trial the channel-index with the maximum z-value
for trllop=1:numtrl
  sgnind = find(maxval(trllop,:) == max(maxval(trllop,:))); 
  artfctchn(trllop) = sgnind;
end;

clear cutfltdatstr;

% update the configuration prior to calling artifact_feedback
cfg.artfctdef.jump = artfctdef;

[artfctdef, artifact] = artifact_feedback_old(cfg, z_tdata, artfctchn, 'jump');

% remember the details that were used here
cfg.artfctdef.jump          = artfctdef;
cfg.artfctdef.jump.artifact = artifact;

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: artifact_jump_old.m,v 1.3 2006/04/20 09:56:07 roboos Exp $';

