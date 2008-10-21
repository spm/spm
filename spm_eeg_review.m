function [] = spm_eeg_review(D,flag)
% function for general review (display) of SPM meeg object
% FORMAT spm_eeg_review(D)
%
% IN:
%   - D: meeg object
%   - flag: (optional) switch to any of the displays.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jean Daunizeau
% $Id: spm_eeg_review.m 2367 2008-10-21 11:00:48Z jean $

D = struct(D);


%-- Initialize SPM figure
D.PSD.handles.hfig = findobj('Tag', 'Graphics');
if isempty(D.PSD.handles.hfig)
    D.PSD.handles.hfig = spm_figure('create','Graphics','Graphics','on');
else
    clf(D.PSD.handles.hfig)
end
set(D.PSD.handles.hfig,'renderer','OpenGL');

%-- Create default userdata structure
[D] = PSD_initUD(D);
if ~strcmp(D.transform.ID,'time')
    D.PSD.type = 'epoched';
    D.PSD.trials.current = 1;
    D.PSD.VIZU.type = 2;
end

%-- Create figure uitabs
labels = {'EEG','MEG','OTHER','info','source'};
callbacks = {'spm_eeg_review_callbacks(''visu'',''main'',''eeg'')',...
    'spm_eeg_review_callbacks(''visu'',''main'',''meg'')',...
    'spm_eeg_review_callbacks(''visu'',''main'',''other'')',...
    'spm_eeg_review_callbacks(''visu'',''main'',''info'')',...
    'spm_eeg_review_callbacks(''visu'',''main'',''source'')'};
try
    [h] = spm_uitab(D.PSD.handles.hfig,labels,callbacks,[],flag);
catch
    clear flag
    [h] = spm_uitab(D.PSD.handles.hfig,labels,callbacks,[],4);
end
D.PSD.handles.tabs = h;

% %-- Initilize display on 'info'
D.PSD.VIZU.modality = 'info';
D.PSD.VIZU.info = 4;
D.PSD.VIZU.fromTab = [];
[D] = spm_eeg_review_switchDisplay(D);

    
%-- Attach userdata to SPM graphics window
set(D.PSD.handles.hfig,'color',[1 1 1],'userdata',D);


if exist('flag','var')
    switch flag
        case 1
            spm_eeg_review_callbacks('visu','main','eeg')
        case 2
            spm_eeg_review_callbacks('visu','main','meg')
        case 3
            spm_eeg_review_callbacks('visu','main','other')
        case 5
            spm_eeg_review_callbacks('visu','main','source')
    end
end




%% initialization of the userdata structure
function [D] = PSD_initUD(D)
% function D = PSD_initUD(D)
% This function initializes the userdata structure.

%-- Check spm_uitable capability (JAVA compatibility) --%
D.PSD.VIZU.uitable = spm_uitable;

%-- Initialize time window basic info --%
D.PSD.VIZU.xlim = [1,min([5e2,D.Nsamples])];
D.PSD.VIZU.info = 1;

%-- Initialize trials info --%
switch D.type
    % before epoching
    case 'continuous'
        D.PSD.type = 'continuous';
        if ~isempty(D.trials) && ~isempty(D.trials(1).events)
            Nevents = length(D.trials.events);
            for i =1:Nevents
                if isempty(D.trials.events(i).duration)
                    D.trials.events(i).duration = 0;
                end
                if isempty(D.trials.events(i).value)
                    D.trials.events(i).value = 0;
                    D.trials.events(i).type = '0';
                end
            end
        end
        D.PSD.VIZU.type = 1;
        
    case 'single' % after epoching
        D.PSD.type = 'epoched';
        nTrials = length(D.trials);
        D.PSD.trials.TrLabels = cell(nTrials,1);
        for i = 1:nTrials
            if D.trials(i).bad
                str = ' (bad)';
            else
                str = ' (not bad)';
            end
            D.PSD.trials.TrLabels{i} = ['Trial ',num2str(i),': ',D.trials(i).label,str];
        end
        D.PSD.trials.current = 1;
        D.PSD.VIZU.type = 1;
    case {'evoked','grandmean'}
        D.PSD.type = 'epoched';
        nTrials = length(D.trials);
        D.PSD.trials.TrLabels = cell(nTrials,1);
        for i = 1:nTrials
            D.PSD.trials.TrLabels{i} = ['Trial ',num2str(i),' (average of ',...
                num2str(D.trials(i).repl),' events): ',D.trials(i).label];
            D.trials(i).events = [];
        end
        D.PSD.trials.current = 1;
        D.PSD.VIZU.type = 1;
end

%-- Initialize channel info --%
nc = length(D.channels);
D.PSD.EEG.I  = find(strcmp('EEG',{D.channels.type}));
D.PSD.MEG.I  = find(strcmp('MEG',{D.channels.type}));
D.PSD.other.I = setdiff(1:nc,[D.PSD.EEG.I(:);D.PSD.MEG.I(:)]);
if ~isempty(D.PSD.EEG.I)
    figure(D.PSD.handles.hfig)
    set(D.PSD.handles.hfig,'userdata',D);
    [out] = spm_eeg_review_callbacks('get','VIZU',D.PSD.EEG.I);
    D.PSD.EEG.VIZU = out;
else
    D.PSD.EEG.VIZU = [];
end
if ~isempty(D.PSD.MEG.I)
    figure(D.PSD.handles.hfig)
    set(D.PSD.handles.hfig,'userdata',D);
    [out] = spm_eeg_review_callbacks('get','VIZU',D.PSD.MEG.I);
    D.PSD.MEG.VIZU = out;
else
    D.PSD.MEG.VIZU = [];
end
if ~isempty(D.PSD.other.I)
    figure(D.PSD.handles.hfig)
    set(D.PSD.handles.hfig,'userdata',D);
    [out] = spm_eeg_review_callbacks('get','VIZU',D.PSD.other.I);
    D.PSD.other.VIZU = out;
else
    D.PSD.other.VIZU = [];
end


%%-- Initialize inverse field info
if isfield(D.other,'inv') && ~isempty(D.other.inv) % && isfield(D.other.inv{1},'inverse')
    isInv = zeros(length(D.other.inv),1);
    for i=1:length(D.other.inv)
        if isfield(D.other.inv{i},'inverse') && strcmp(D.other.inv{i}.method,'Imaging')
            isInv(i) = 1;
        end
    end
    isInv = find(isInv);
    Ninv = length(isInv);
    if Ninv>=1
        labels = cell(Ninv,1);
        callbacks = cell(Ninv,1);
        F = zeros(Ninv,1);
        pst = [];
        for i=1:Ninv
            if ~isfield(D.other.inv{isInv(i)},'comment')
                D.other.inv{isInv(i)}.comment{1} = num2str(i);
            end
            if ~isfield(D.other.inv{isInv(i)},'date')
                D.other.inv{isInv(i)}.date(1,:) = 'unknown';
                D.other.inv{isInv(i)}.date(2,:) = '       ';
            end
            labels{i} = [D.other.inv{isInv(i)}.comment{1}];
            callbacks{i} = ['spm_eeg_review_callbacks(''visu'',''inv'',',num2str(i),')'];
            F(i) = D.other.inv{isInv(i)}.inverse.F;
            pst = [pst;D.other.inv{isInv(i)}.inverse.pst(:)];
        end
        pst = unique(pst);
    end
else
    Ninv = 0;
end
if Ninv >= 1
    D.PSD.source.VIZU.current = 1;
    D.PSD.source.VIZU.isInv = isInv;
    D.PSD.source.VIZU.pst = pst;
    D.PSD.source.VIZU.F = F;
    D.PSD.source.VIZU.labels = labels;
    D.PSD.source.VIZU.callbacks = callbacks;
    D.PSD.source.VIZU.timeCourses = 1;
else
    D.PSD.source.VIZU.current = 0;
    D.PSD.source.VIZU.isInv = [];
    D.PSD.source.VIZU.pst = [];
    D.PSD.source.VIZU.F = [];
    D.PSD.source.VIZU.labels = [];
    D.PSD.source.VIZU.callbacks = [];
    D.PSD.source.VIZU.timeCourses = [];
end




