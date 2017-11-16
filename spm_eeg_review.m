function spm_eeg_review(D,flag,inv)
% General review (display) of SPM meeg object
% FORMAT spm_eeg_review(D,flags,inv)
%
% INPUT:
% D      - meeg object
% flag   - switch to any of the displays (optional)
% inv    - which source reconstruction to display (when called from
% spm_eeg_inv_imag_api.m)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jean Daunizeau
% $Id: spm_eeg_review.m 7221 2017-11-16 14:25:37Z vladimir $

if nargin == 0
    [D, sts] = spm_select(1, 'mat$', 'Select M/EEG mat file');
    if ~sts, return; end
    D = spm_eeg_load(D);
end

%-- Initialize SPM figure
D.PSD.handles.hfig = spm_figure('GetWin','Graphics');
spm_clf(D.PSD.handles.hfig);
% Get default SPM graphics options --> revert back to defaults
D.PSD.SPMdefaults.col = get(D.PSD.handles.hfig,'colormap');
D.PSD.SPMdefaults.renderer = get(D.PSD.handles.hfig,'renderer');

%-- Create default userdata structure
try
    D.PSD.source.VIZU.current = inv;
end
[D] = PSD_initUD(D);
if ~strcmp(transformtype(D),'time') % 
    D.PSD.type = 'epoched';
    D.PSD.trials.current = 1;
    D.PSD.VIZU.type = 2;
end

%-- Create figure uitabs
labels = {'EEG', 'MEG', 'MPLANAR', 'MCOMB', 'OTHER','info','source'};
callbacks = {'spm_eeg_review_callbacks(''visu'',''main'',''eeg'')',...
    'spm_eeg_review_callbacks(''visu'',''main'',''meg'')',...
    'spm_eeg_review_callbacks(''visu'',''main'',''megplanar'')',...
    'spm_eeg_review_callbacks(''visu'',''main'',''megcomb'')',...
    'spm_eeg_review_callbacks(''visu'',''main'',''other'')',...
    'spm_eeg_review_callbacks(''visu'',''main'',''info'')',...
    'spm_eeg_review_callbacks(''visu'',''main'',''source'')'};
try
    [h] = spm_uitab(D.PSD.handles.hfig,labels,callbacks,[],flag);
catch
    [h] = spm_uitab(D.PSD.handles.hfig,labels,callbacks,[],6);
end
D.PSD.handles.tabs = h;

% Add prepare and SAVE buttons
object.type = 'buttons';
object.list = 1;
D = spm_eeg_review_uis(D,object);
set(D.PSD.handles.BUTTONS.pop1,...
    'deletefcn',@back2defaults)

%-- Attach userdata to SPM graphics window
Dtmp = rmfield(D,'PSD');
D.PSD.D0 = Dtmp;
%CP, some how, removing and adding a field in the object doesn't work directly...
set(D.PSD.handles.hfig,...
    'units','normalized',...
    'color',[1 1 1],...
    'userdata',D);

try
    if ismac
        set(D.PSD.handles.hfig,'renderer','zbuffer');
    else
        set(D.PSD.handles.hfig,'renderer','OpenGL');
    end
catch
    set(D.PSD.handles.hfig,'renderer','OpenGL');
end

try
    switch flag
        case 1
            spm_eeg_review_callbacks('visu','main','eeg')
        case 2
            spm_eeg_review_callbacks('visu','main','meg')
        case 3
            spm_eeg_review_callbacks('visu','main','megplanar')
        case 4
            spm_eeg_review_callbacks('visu','main','megcomb')
        case 5
            spm_eeg_review_callbacks('visu','main','other')
        case 6
            spm_eeg_review_callbacks('visu','main','info')
        case 7
            spm_eeg_review_callbacks('visu','main','source')
    end
catch
    % Initilize display on 'info'
    spm_eeg_review_callbacks('visu','main','info')
end


%% Revert graphical properties of SPM figure back to normal
function back2defaults(e1,e2)
hf = spm_figure('FindWin','Graphics');
D = get(hf,'userdata');
try
    set(D.PSD.handles.hfig,'colormap',D.PSD.SPMdefaults.col);
    set(D.PSD.handles.hfig,'renderer',D.PSD.SPMdefaults.renderer);
end


%% initialization of the userdata structure
function [D] = PSD_initUD(D)
% This function initializes the userdata structure.

%-- Check spm_uitable capability (JAVA compatibility) --%
D.PSD.VIZU.uitable = spm_uitable;

%-- Initialize time window basic info --%
D.PSD.VIZU.xlim = [1,min([5e2,nsamples(D)])];
D.PSD.VIZU.info = 4; % show history
D.PSD.VIZU.fromTab = [];


%-- Initialize trials info --%
switch D.type
    
    %------ before epoching -----%
    
    case 'continuous'
        D.PSD.type = 'continuous';
        if ntrials(D) && ~isempty(events(D,1))
            Events = events(D);
            Nevents = length(Events);
            for i =1:Nevents
                if isempty(Events(i).duration)
                    Events(i).duration = 0;
                end
                if isempty(Events(i).value)
                    Events(i).value = '0';
                end
                if isempty(Events(i).type)
                    Events(i).type = '0';
                end
                if ~ischar(Events(i).value)
                    Events(i).value = num2str(Events(i).value);
                end
                if ~ischar(Events(i).type)
                    Events(i).type = num2str(Events(i).type);
                end
            end
            %D = events(D,1,Events);
        end
        D.PSD.VIZU.type = 1;
        
        %------ after epoching -----%
        
    case 'single'
        D.PSD.type = 'epoched';
        nTrials = D.ntrials;
        D.PSD.trials.TrLabels = cell(nTrials,1);
        bTrials = badtrials(D);
        for i = 1:nTrials
            if any(i==bTrials)
                str = ' (bad)';
            else
                str = ' (not bad)';
            end
            D.PSD.trials.TrLabels{i} = [...
                'Trial ',num2str(i),': ',char(conditions(D,i)) ,str];
        end
        D.PSD.trials.current = 1;
        D.PSD.VIZU.type = 1;
    case {'evoked','grandmean'}
        D.PSD.type = 'epoched';
        nTrials = D.ntrials;
        D.PSD.trials.TrLabels = cell(nTrials,1);
        for i = 1:nTrials
            D.PSD.trials.TrLabels{i} = [...
                'Trial ',num2str(i),' (average of ',...
                num2str(repl(D,i)),' events): ',...
                char(conditions(D,i))];
            D = events(D,i,[]);
        end
        D.PSD.trials.current = 1;
        D.PSD.VIZU.type = 1;
end

%-- Initialize channel info --%
nc = D.nchannels;
D.PSD.EEG.I  = indchantype(D,'EEG');
D.PSD.MEG.I  = sort(indchantype(D,'MEG'));
D.PSD.MEGPLANAR.I  = indchantype(D,'MEGPLANAR');
D.PSD.MEGCOMB.I  = indchantype(D,'MEGCOMB');
D.PSD.other.I = setdiff(1:nc, ...
    [D.PSD.EEG.I(:);D.PSD.MEG.I(:);D.PSD.MEGPLANAR.I(:);D.PSD.MEGCOMB.I(:)]);
%-- Get basic display variables (data range, offset,...)
if ~isempty(D.PSD.EEG.I)
    set(D.PSD.handles.hfig,'userdata',D);
    figure(D.PSD.handles.hfig)
    [out] = spm_eeg_review_callbacks('get','VIZU',D.PSD.EEG.I);
    D.PSD.EEG.VIZU = out;
else
    D.PSD.EEG.VIZU = [];
end
if ~isempty(D.PSD.MEG.I)
    set(D.PSD.handles.hfig,'userdata',D);
    figure(D.PSD.handles.hfig)
    [out] = spm_eeg_review_callbacks('get','VIZU',D.PSD.MEG.I);
    D.PSD.MEG.VIZU = out;
else
    D.PSD.MEG.VIZU = [];
end
if ~isempty(D.PSD.MEGPLANAR.I)
    set(D.PSD.handles.hfig,'userdata',D);
    figure(D.PSD.handles.hfig)
    [out] = spm_eeg_review_callbacks('get','VIZU',D.PSD.MEGPLANAR.I);
    D.PSD.MEGPLANAR.VIZU = out;
else
    D.PSD.MEGPLANAR.VIZU = [];
end
if ~isempty(D.PSD.MEGCOMB.I)
    set(D.PSD.handles.hfig,'userdata',D);
    figure(D.PSD.handles.hfig)
    [out] = spm_eeg_review_callbacks('get','VIZU',D.PSD.MEGCOMB.I);
    D.PSD.MEGCOMB.VIZU = out;
else
    D.PSD.MEGCOMB.VIZU = [];
end
if ~isempty(D.PSD.other.I)
    set(D.PSD.handles.hfig,'userdata',D);
    figure(D.PSD.handles.hfig)
    [out] = spm_eeg_review_callbacks('get','VIZU',D.PSD.other.I);
    D.PSD.other.VIZU = out;
else
    D.PSD.other.VIZU = [];
end


%-- Initialize inverse field info --%
if isfield(D,'inv') && ~isempty(D.inv)
    isInv = zeros(length(D.inv),1);
    for i=1:length(D.inv)
        if isfield(D.inv{i},'inverse') && ...
                isfield(D.inv{i}, 'method') && ...
                strcmp(D.inv{i}.method,'Imaging')
            isInv(i) = 1;
        end
    end
    isInv = find(isInv);
    Ninv = length(isInv);
    if Ninv>=1
        labels = cell(Ninv,1);
        callbacks = cell(Ninv,1);
        F  = zeros(Ninv,1);
        ID = zeros(Ninv,1);
        pst = [];
        for i=1:Ninv
            if ~isfield(D.inv{isInv(i)},'comment')
                D.inv{isInv(i)}.comment{1} = num2str(i);
            end
            if ~isfield(D.inv{isInv(i)},'date')
                D.inv{isInv(i)}.date(1,:) = '?';
                D.inv{isInv(i)}.date(2,:) = ' ';
            end
            if isfield(D.inv{isInv(i)}.inverse,'R2') ...
                    && isnan(D.inv{isInv(i)}.inverse.R2)
                D.inv{isInv(i)}.inverse.R2 = [];
            end
            if isfield(D.inv{isInv(i)}.inverse, 'ID')
                ID(i) = D.inv{isInv(i)}.inverse.ID;
            else
                ID(i) = nan;
            end
            labels{i} = [D.inv{isInv(i)}.comment{1}];
            callbacks{i} = ['spm_eeg_review_callbacks(''visu'',''inv'',',num2str(i),')'];
            try
                F(i) = D.inv{isInv(i)}.inverse.F;
                pst = [pst;D.inv{isInv(i)}.inverse.pst(:)];
            catch
                continue
            end
        end
        if isempty(pst)
            Ninv = 0;
        else
            pst = unique(pst);
        end
    end
else
    Ninv = 0;
end
if Ninv >= 1
    try
        if D.PSD.source.VIZU.current > Ninv
            D.PSD.source.VIZU.current = 1;
        end
    catch
        D.PSD.source.VIZU.current = 1;
    end
    D.PSD.source.VIZU.isInv = isInv;
    D.PSD.source.VIZU.pst = pst;
    D.PSD.source.VIZU.F = F;
    D.PSD.source.VIZU.ID = ID;
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
