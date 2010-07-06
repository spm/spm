function [stats,talpositions]=spm_eeg_ft_beamformer_gui(S)
% LCMV univariate beamformer
% FORMAT [stats,talpositions]=spm_eeg_ft_beamformer_gui(S)
% 
% S            - struct (optional)
% (optional) fields of S:
% S.D          - meeg object or filename
%                coregistration must has been performed beforehand
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Gareth Barnes
% $Id: spm_eeg_ft_beamformer_gui.m 3971 2010-07-06 09:53:38Z gareth $

[Finter,Fgraph] = spm('FnUIsetup','LCMV beamformer for power', 0);
%%

%% ============ Load SPM EEG file and verify consistency
if nargin == 0
    S = [];
end

try
    D = S.D;
catch
    D = spm_select(1, '\.mat$', 'Select MEEG mat file');
    S.D = D;
end

if ischar(D)
    try
        D = spm_eeg_load(D);
    catch
        error(sprintf('Trouble reading file %s', D));
    end
end

[ok, D] = check(D, 'sensfid');

if ~ok
    if check(D, 'basic')
        errordlg(['The requested file is not ready for source reconstruction.'...
            'Use prep to specify sensors and fiducials.']);
    else
        errordlg('The meeg file is corrupt or incomplete');
    end
    return
end

modality = spm_eeg_modality_ui(D, 1, 1);
channel = D.chanlabels(strmatch(modality, D.chantype))';

if isfield(S, 'refchan') && ~isempty(S.refchan)
    refchan = S.refchan;
else
    refchan = [];
end

%% ============ Find or prepare head model

if ~isfield(D, 'val')
    D.val = 1;
end



[ok, D] = check(D, 'sensfid');

if ~ok
    if check(D, 'basic')
        errordlg(['The requested file is not ready for source reconstruction.'...
            'Use prep to specify sensors and fiducials.']);
    else
        errordlg('The meeg file is corrupt or incomplete');
    end
    return
end


 if ~isfield(D,'inv')
     errordlg('Need to set up a forward model before you start');
     return;
 end;



clb = D.condlist;
triallist=[];
trialtypes=[];
latencies=[];
outfilenames='';

contrast_str='';
if ~isfield(S,'timewindows'),
    S.timewindows={[],[]};
end; % if

Ntimewindows=2; % %% fixed for now
for i = 1:Ntimewindows, %%spm_input('Number of time windows:', '+1', 'r', '1', 1)
    if i==1,
        activestr='Active condition (+1)'
    else
        activestr='Baseline condition (-1)'
    end;
    contrast_str=[contrast_str 'Tr']
    if isempty(S.timewindows{i});
        [selection, ok]= listdlg('ListString', clb, 'SelectionMode', 'multiple' ,'Name', sprintf('Select %s',activestr) , 'ListSize', [400 300]);
        if ~ok
            return;
        end
        S.trigger{i}=clb(selection);
        outstr=sprintf('Offset (ms) from %s',cell2mat(S.trigger{i}));
        starttime= spm_input(outstr, '+1', 'r', '', 1);
        if i==1,
            outstr=sprintf('Duration (ms) ');
            duration= spm_input(outstr, '+1', 'r', '', 1);
        end; % if i==1
        S.timewindows{i}=[starttime;starttime+duration]/1000;
    else
        duration=(S.timewindows{1}(2)-S.timewindows{1}(1))*1000; %% convert back to ms
        duration_check=(S.timewindows{2}(2)-S.timewindows{2}(1))*1000;
        if duration~=duration_check,
            error('durations of time windows not equal');
        end;
        
    end; % if isfield timewindows
    
    Ntrialspercond=0;
    
    for s=1:numel(S.trigger{i}), %% maybe multiple triggers in selection
        trigname=cell2mat(S.trigger{i}(s));
        Ntrialspercond=Ntrialspercond+length(D.pickconditions(trigname)); %% number of trials for this condition
        contrast_str=[contrast_str trigname];
        triallist=[triallist , D.pickconditions(trigname)];
        contrast_str=[contrast_str ','];
    end; % for selection
    trialtypes=[trialtypes ; ones(Ntrialspercond,1).*i];
    latencies =[latencies;repmat([S.timewindows{i}],1,Ntrialspercond)'];
    contrast_str=[contrast_str num2str(S.timewindows{i}(1)) '-' num2str(S.timewindows{i}(2)) 's'];
    if i==1,
        contrast_str=[contrast_str '_vs'];
    end; % if i
end; % for i

contrast_str=sprintf('%s',deblank(contrast_str));

Nconditions=numel(S.timewindows);

type1ind=find(trialtypes==1);
type2ind=find(trialtypes==2);

if 2*abs((length(type1ind)-length(type2ind)))./(length(type1ind)+length(type2ind))>0.1,

    balance = spm_input('trial numbers diffe (>10%). Randomly resample ?','+1', 'yes|no', [1, 0]);
    minlen=min(length(type1ind),length(type2ind));
    m1=randperm(length(type1ind));
    type1ind=sort(type1ind(m1(1:minlen)));
    m1=randperm(length(type2ind));
    type2ind=sort(type2ind(m1(1:minlen)));
    disp(sprintf('Now both conditions have %d trials',minlen));
end;

    
    

%% Set up design matrix for a t test
S.design.X=size(latencies,1);
S.design.X(type1ind,1)=1;
S.design.X(type2ind,2)=1;
S.design.X(:,3)=1;


 
 contrast=[1 -1 0];
      S.design.contrast=contrast;
      S.design.Xwindowduration=duration/1000; %% in seconds 
      S.design.Xtrials=triallist'; % correspond to the trials 
      S.design.Xstartlatencies=latencies(:,1); %% correspond to start latencies within trials


freqbands=[];
if ~isfield(S, 'freqbands')
    for i = 1:spm_input('Number of frequency bands:', '+1', 'r', '1', 1)
        outstr=sprintf('Band %d [st end] in Hz ',i);
        S.freqbands{i} = spm_input(outstr, '+1', 'r', '', 2)';
        freqbands =[freqbands;S.freqbands{i}];
    end
end;

Nbands=numel(S.freqbands);

if ~isfield(S,'gridstep');
 S.gridstep = spm_input('Grid step (mm):', '+1', 'r', '5');
end; 


if ~isfield(S,'logflag'),
    S.logflag=[];
end; % if

if isempty(S.logflag),
    S.logflag=0;
end; % if
  
 if ~isfield(S, 'regpc')
     S.regpc =spm_input('Regularization (%):', '+1', 'r', '0');
 end

if ~isfield(S, 'preview')
    S.preview = spm_input('Preview results?','+1', 'yes|no', [1, 0]);
end
%%
 
 [stats,talpositions]=spm_eeg_ft_beamformer_lcmv(S);
 
    
end % function
