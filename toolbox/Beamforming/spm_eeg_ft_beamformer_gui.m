function [stats,talpositions]=spm_eeg_ft_beamformer_gui(S)
% _______________________________________________________________________
% Copyright (C) 2009 Institute of Neurology, UCL
% basic gui for an LCMV univariate beamformer
%
% Gareth Barnes
% $Id: spm_eeg_ft_beamformer_gui.m 3514 2009-10-28 14:37:09Z gareth $

[Finter,Fgraph] = spm('FnUIsetup','LCMV beamformer for power', 0);
%%

%% ============ Load SPM EEG file and verify consistency
if nargin == 0
    S = [];
end

try
    D = S.D;
catch
    D = spm_select(1, '\.mat$', 'Select EEG mat file');
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



for m = 1:numel(D.inv{D.val}.forward)
    if strncmp(modality, D.inv{D.val}.forward(m).modality, 3)
        vol  = D.inv{D.val}.forward(m).vol;
        if isa(vol, 'char')
            vol = fileio_read_vol(vol);
        end
        datareg  = D.inv{D.val}.datareg(m);
    end
end

 try
     vol = D.inv{D.val}.forward.vol;
     datareg = D.inv{D.val}.datareg;
 catch
     D = spm_eeg_inv_mesh_ui(D, D.val, [], 1);
     D = spm_eeg_inv_datareg_ui(D, D.val);
     datareg = D.inv{D.val}.datareg;
 end


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


%% Set up design matrix for a t test
S.design.X=size(latencies,1);
S.design.X(find(trialtypes==1),1)=1;
S.design.X(find(trialtypes==2),1)=-1;
 
 contrast=[1];
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


