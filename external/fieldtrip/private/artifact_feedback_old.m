function [artfctdef, artifact] = artifact_feedback(cfg,z_tdata,artfctchn,flag);

% ARTIFACT_FEEDBACK provides feedback about the artifacts in the data.
%
% [artfctdef, artifact] = artifact_feedback(cfg,artfctdef,artifact,z_tdata,artfctchn);
% If cfg.artfctdef.xxx.feedback was set to 'yes', the function updates the artifact-definition
% parameters and outputs artifact. artifact is a binary vector as a function of time (in samples).
% Its value is 1 if in one of the channels to be analysed an xxx artifact
% is present. artifact is computed by thresholding the vector z_tdata, and
% convolving this result with an appropriate artifact-kernel.
% The rejection threshold can be set interactively, while given feedback about the to
% be rejected epochs. This can be done, until the threshold is set in a satisfactory way.

% Copyright (c) 2003,  F.C. Donders Centre
%
% $Log: artifact_feedback_old.m,v $
% Revision 1.1  2006/01/12 14:21:26  roboos
% renamed JMs artifact_feedback and artifact_viewer into xxx_old
% they are only provided for people who insist on running the old code
%
% Revision 1.14  2005/08/01 12:40:22  roboos
% added units (samples/uV) to the axes in the feedback plots
%
% Revision 1.13  2005/05/17 17:50:49  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.12  2005/04/22 08:00:14  jansch
% fixed bug that caused routine to crash in the case of no artifacts
%
% Revision 1.11  2005/04/19 16:50:36  roboos
% do not try to prevent artifact leakage if there are no artifacts
%
% Revision 1.10  2005/04/19 12:33:51  jansch
% fixed minor bug
%
% Revision 1.9  2005/04/08 12:07:09  roboos
% added an empty line between the help and the log messages
%
% Revision 1.8  2005/04/08 11:54:30  jansch
% notwithstanding pieter's excellent feedback, Robert did not manage to fix
% the bug sufficiently. now it should work
%
% Revision 1.7  2005/04/08 11:47:23  jansch
% fixed small bug.
%
% Revision 1.6  2005/04/08 11:33:13  roboos
% prolonged the rejectall vector, so that it can hold all trials after artifact padding (thanks to Pieter Medendorp)
%
% Revision 1.5  2005/04/08 08:23:47  roboos
% replaced convolve_suprathreshold subfunction with two separate functions, one that does the padding (not using a convolution) and that merges the artifacts, the other one preventing the artifacts from leaking into other trials for non-continuous data
%
% Revision 1.4  2005/04/08 07:10:09  roboos
% added a check to prevent a crash due to a bug (negative artifact length), the bug is still present
%
% Revision 1.3  2005/04/06 13:08:00  roboos
% implemented check for non-continuous data to ensure that artifacts do not leak over trial boundaries
%
% Revision 1.2  2005/02/08 08:36:55  jansch
% fixed bug in handling of non-padded data
%
% Revision 1.1  2005/02/07 15:10:50  roboos
% moved two subfunctions for artifact detection from main to private directory
%
% Revision 1.4  2004/09/01 09:39:28  jansch
% *** empty log message ***
%
% Revision 1.3  2004/09/01 09:27:25  jansch
% fixed some small bugs, and inconsistencies. replaced convolution with fft-convolution in the merging of artifacts, closely related in time. added GUI in artifact_viewer
%

if     strcmp(flag,'muscle'),    artfctdef = cfg.artfctdef.muscle; 
elseif strcmp(flag,'jump'),      artfctdef = cfg.artfctdef.jump;
elseif strcmp(flag,'EOGraw') || strcmp(flag,'EOGflt') || strcmp(flag,'EOG') ,artfctdef = cfg.artfctdef.eog;
end;

hdr       = read_fcdc_header(cfg.headerfile);
% artifact  = convolve_suprathreshold(cfg,artfctdef,z_tdata,hdr.Fs);
artifact  = artifact_padding(cfg,artfctdef,z_tdata,hdr.Fs);
artifact  = prevent_trial_leakage(cfg,artfctdef,artifact,hdr.Fs);


%rejectall = zeros(1,max(numtrlsmp, numartsmp));	% make a vector with enough samples to hold the trials and artifacts
rejectall = zeros(1,length(z_tdata));
if ~isempty(artifact),
  % update the rejectall vector
  numtrlsmp = max(cfg.trl(:,2));			% count the number of samples in the trial definition
  numartsmp = max(artifact(:,2));			% count the number of samples in the trial definition, could be longer since the artifact can be padded
  for i=1:size(artifact,1)
    rejectall(artifact(i,1):artifact(i,2)) = 1;
  end
end

if strcmp(artfctdef.feedback,'yes') || strcmp(artfctdef.feedback,'flt')
	
	trl    = cfg.trl;
	numtrl = size(trl,1);
	accept = 0; 
	
	while accept == 0 ,         
           h = figure;
           plot(z_tdata);
           hold on;
           plot([1 length(z_tdata)],[artfctdef.cutoff artfctdef.cutoff],'r:');
           hold off;
           xlabel('samples');
           ylabel('zscore');
   	
           response = input('\nwould you like to page through the data (y/n) ?\n','s');
           switch response
           case 'y'
               close(h);
               artifact_viewer_old(cfg,artfctdef,rejectall,z_tdata,artfctchn,flag);
           end;
   	
           fprintf(strcat('\ncurrent  ',artfctdef.method,' threshold = %1.3f'),artfctdef.cutoff);
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
           if oldcutoff ~= artfctdef.cutoff,
   	      % artifact = convolve_suprathreshold(cfg,artfctdef,z_tdata,hdr.Fs);
              artifact  = artifact_padding(cfg,artfctdef,z_tdata,hdr.Fs);
              artifact  = prevent_trial_leakage(cfg,artfctdef,artifact,hdr.Fs);
              % update the rejectall vector
              rejectall(:) = 0;
              for i=1:size(artifact,1)
                rejectall(artifact(i,1):artifact(i,2)) = 1;
              end
           end
           close;
	end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [artifact] = artifact_padding(cfg,artfctdef,z_tdata,Fs);
fprintf('\nidentifying artifact-epochs\n');
datind = z_tdata>artfctdef.cutoff;
begsmp = (find(diff([0 datind])==1))';
endsmp = (find(diff([datind 0])==-1))';
% add the padding to the begin and end of each artifact
begsmp = begsmp - round(artfctdef.pretim .* Fs);
endsmp = endsmp + round(artfctdef.psttim .* Fs);
begsmp(begsmp<1) = 1;       % prevent artifacts to start before time zero
% merge the artifacts that are overlapping each other
remove = zeros(size(begsmp));
for i=2:length(begsmp)
  if begsmp(i)<=endsmp(i-1)
    begsmp(i)= begsmp(i-1);	% include the previous artifact into this one
    remove(i-1) = 1;		% mark the previous one so that it will be removed
  end
end
begsmp = begsmp(~remove);	% remove the artifacts that have been swallowed up
endsmp = endsmp(~remove);	% remove the artifacts that have been swallowed up
% assign the output artifact structure
if ~isempty(begsmp) && ~isempty(endsmp),
  artifact = [begsmp(:) endsmp(:)];
else
  artifact = [];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [artifact] = prevent_trial_leakage(cfg,artfctdef,artifact,Fs);
if ~strcmp(cfg.datatype, 'continuous')
  if isempty(artifact)
    fprintf('no artifacts are found, therefore no leakage will occur\n');
    return
  end
  fprintf('preventing artifacts from leaking into other trials\n');
  trl = cfg.trl;
  dum = artifact;
  % remove the padding at the begin and end of each artifact
  dum(:,1) = dum(:,1) + round(artfctdef.pretim*Fs);
  dum(:,2) = dum(:,2) - round(artfctdef.psttim*Fs);
  for i=1:size(dum,1);
   containcomplete = trl(:,1)<=dum(i,1) & trl(:,2)>=dum(i,2);
   % containpartial  = trl(:,1)<=dum(i,2) & trl(:,2)>=dum(i,1);
   % containpartial  = containpartial & ~containcomplete;
   sel = find(containcomplete);
   if length(sel)==1
      artifact(i,1) = max(trl(sel,1), artifact(i,1));
      artifact(i,2) = min(trl(sel,2), artifact(i,2));
   elseif length(sel)>1
     error('artifact is completely contained within two trials, error in trial definition');
   end
  end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [artifact, datind] = convolve_suprathreshold(cfg,artfctdef,z_tdata,Fs);

fprintf('\nidentifying artifact-epochs\n');

trl    = cfg.trl;
presmp = max([round(artfctdef.pretim .* Fs) 1]);
pstsmp = round(artfctdef.psttim .* Fs);
maxsmp = max([presmp,pstsmp]);
kernel = [zeros(1,maxsmp - presmp) ones(1,presmp) ...
          ones(1,pstsmp) zeros(1,maxsmp - pstsmp)];
krnsmp = length(kernel);
datsmp = size(z_tdata,2);
datind = zeros(size(z_tdata));
datind = sum(z_tdata>artfctdef.cutoff,1);
%datind = (conv(datind,kernel))>0;

begsmp = (find(diff([0 datind])==1))';
endsmp = (find(diff([datind 0])==-1))';

% add the padding to the begin and end of each artifact
begsmp = begsmp - presmp;
endsmp = endsmp - pstsmp;
begsmp(begsmp<1) = 1;       % prevent artifacts before time zero

% merge overlapping artifacts together
remove = zeros(size(begsmp));
for i=2:length(begsmp)
  if begsmp(i)<=endsmp(i-1)
    begsmp(i)= begsmp(i-1);	% merge with the previous one
    remove(i-1) = 1;		% mark the previous one so that it will be removed
  end
end
begsmp = begsmp(~remove);	% remove the artifacts that have been swallowed up
endsmp = endsmp(~remove);	% remove the artifacts that have been swallowed up

% assign the output artifact structure
if ~isempty(begsmp) && ~isempty(endsmp),
    artifact(:,1) = begsmp;
    artifact(:,2) = endsmp;
else
    artifact = [];
end;

if ~strcmp(cfg.datatype, 'continuous')
  fprintf('preventing artifacts from leaking into other trials\n');
  dum = artifact;
  dum(:,1) = dum(:,1) + round(artfctdef.pretim*Fs);
  dum(:,2) = dum(:,2) - round(artfctdef.psttim*Fs);
  % somehow the convolution does not do what I would expect
  % the pretim and psttim are not added correctly around the artifact
  sel = dum(:,1)>dum(:,2);	% find artifacts with negative length
  dum(sel,2) = dum(sel,1);	% set the length to zero
  for i=1:size(dum,1);
   containcomplete = trl(:,1)<=dum(i,1) & trl(:,2)>=dum(i,2);
   % containpartial  = trl(:,1)<=dum(i,2) & trl(:,2)>=dum(i,1);
   % containpartial  = containpartial & ~containcomplete;
   sel = find(containcomplete);
   if length(sel)==1
      artifact(i,1) = max(trl(sel,1), artifact(i,1));
      artifact(i,2) = min(trl(sel,2), artifact(i,2));
   elseif length(sel)>1
     error('this is weird, contact Robert');
   end
  end
end

