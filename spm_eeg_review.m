function [] = spm_eeg_review(D,flag)
% function for general review (display) of SPM meeg object
% FORMAT spm_eeg_review(D)
%
% IN:
%   - D: meeg object
%   - flag: (optional) switch to any of the displaying ('standardData',
%   'scalpData' or 'visuRecon' -> default is 'standardData').
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jean Daunizeau
% $Id: spm_eeg_review.m 1826 2008-06-16 13:51:36Z guillaume $

if ~exist('flag','var') || isempty(flag)
    flag = 'standardData';
    active = 1;
else
    switch flag
        case 'standardData'
            active = 1;
        case 'scalpData'
            active = 2;
        case 'visuRecon'
            active = 3;
        otherwise
            flag = 'standardData';
            active = 1;
    end

end

%----------- Update MATLAB path -------------------------%
% [PSD_rootPath] = fileparts(mfilename('fullpath'));
% addpath(PSD_rootPath)
% addpath([PSD_rootPath,filesep,'routinetheque'])
% addpath([PSD_rootPath,filesep,'templates'])
% addpath([PSD_rootPath,filesep,'openFilters'])

%----------- Create default userdata structure ----------%
[D]                         = PSD_initUD(D);


%----------Initialize display stuff -----------------------%
D.PSD.handles.hfig          = findobj('Tag', 'Graphics');
if isempty(D.PSD.handles.hfig)
    D.PSD.handles.hfig      = spm_figure('create','Graphics','Graphics','on');
else
    clf(D.PSD.handles.hfig)
end
set(D.PSD.handles.hfig,'renderer','OpenGL');
% POS = get(D.PSD.handles.hfig,'position');
% D.PSD.handles.hfig = figure('position',POS);
% D.PSD.handles.Rotate3d = rotate3d;
% set(D.PSD.handles.Rotate3d,'enable','off');
% clf(D.PSD.handles.hfig)
% set(D.PSD.handles.hfig, 'SelectionType', 'normal');


%---------- Create GUI usermenu and get the handles ----------%
[D.PSD.handles]             = createPSDuimenu(D,1);
% switch D.PSD.type
%     case 'epoched'
%         labels = {'standard','scalp','source reconstruction'};
%         callbacks = {'spm_eeg_review_callbacks(''visu'',''main'',1)',...
%             'spm_eeg_review_callbacks(''visu'',''main'',2)',...
%             'spm_eeg_review_callbacks(''visu'',''main'',3)'};
%     case 'continuous'
%         labels = {'standard'};
%         callbacks = [];
% end
labels = {'standard','scalp','source'};
callbacks = {'spm_eeg_review_callbacks(''visu'',''main'',1)',...
    'spm_eeg_review_callbacks(''visu'',''main'',2)',...
    'spm_eeg_review_callbacks(''visu'',''main'',3)'};
[h] = spm_uitab(D.PSD.handles.hfig,labels,callbacks,[],active);
D.PSD.handles.tabs = h;
D.PSD.VIZU.type             = [];
[D]                         = spm_eeg_review_switchDisplay(D,flag);


set(D.PSD.handles.hfig,'color',[1 1 1]);
set(D.PSD.handles.hfig,'userdata',D);

if ~isempty(D.PSD.VIZU.type)
    spm_eeg_review_callbacks('visu','update')
    if strcmp(D.PSD.VIZU.type,'standardData')
        spm_eeg_review_callbacks('visu','time_w',0.5)
    end
else
    D.PSD.VIZU.type = flag;
end





%% Create ui menus and buttons
function [handles] = createPSDuimenu(D,flag)
% This function creates the UI menu associated with the preSelectData toolbox.
% IN:
%   - D: userdata structure
%   - flag: binary variable (=1 if data is provided, =0 otherwise)
% OUT:
%   - handles: structure containing the handles of all controls created for
%   the preSelectData GUI.

if flag == 1
    enab = 'on';
else
    enab = 'off';
end

hfig = D.PSD.handles.hfig;
handles = D.PSD.handles;
figure(hfig)
% POS = get(hfig,'position');
%
% switch D.PSD.type
%     case 'epoched'
%         labels = {'standard','scalp','source reconstruction'};
%         callbacks = {'spm_eeg_review_callbacks(''visu'',''main'',1)',...
%     'spm_eeg_review_callbacks(''visu'',''main'',2)',...
%     'spm_eeg_review_callbacks(''visu'',''main'',3)'};
%     case 'continuous'
%         labels = {'standard'};
%         callbacks = [];
% end
% [h] = spm_uitab(hfig,labels,callbacks);
% handles.tabs = h;

% % VISUALIZATION uimenu
% f = uimenu('Label','Visualization','enable',enab);
% handles.VIZU.scalp = uimenu(f,'Label','Scalp interpolation');
% handles.VIZU.scalp.interp1 = uimenu(handles.VIZU.scalp,'Label','Image scalp data','Callback','spm_eeg_review_callbacks(''visu'',''scalp_interp'',1)');
% handles.VIZU.zoom = uimenu(f,'Label','Zoom');
% handles.VIZU.zoom1 = uimenu(handles.VIZU.zoom,'Label','Zoom in (mouse box)','Callback','spm_eeg_review_callbacks(''visu'',''zoom'',1)');
% handles.VIZU.zoom2 = uimenu(handles.VIZU.zoom,'Label','Zoom reset','Callback','spm_eeg_review_callbacks(''visu'',''zoom'',0)','enable','off');
% handles.VIZU.inten_sc = uimenu(f,'Label','Intensity rescaling');
% handles.VIZU.inten_sc1 = uimenu(handles.VIZU.inten_sc,'Label','Contrast +++','Callback','spm_eeg_review_callbacks(''visu'',''iten_sc'',2)');
% handles.VIZU.inten_sc2 = uimenu(handles.VIZU.inten_sc,'Label','Contrast ---','Callback','spm_eeg_review_callbacks(''visu'',''iten_sc'',0.5)');
% handles.VIZU.time_w = uimenu(f,'Label','Size of the plotted time window');
% handles.VIZU.time_w1 = uimenu(handles.VIZU.time_w,'Label','Width +++','Callback','spm_eeg_review_callbacks(''visu'',''time_w'',2)');
% handles.VIZU.time_w2 = uimenu(handles.VIZU.time_w,'Label','Width ---','Callback','spm_eeg_review_callbacks(''visu'',''time_w'',0.5)');
% handles.VIZU.f_sensors = uimenu(f,'Label','Choose channels');
% handles.VIZU.visu_sensor_select = uimenu(handles.VIZU.f_sensors,'Label','Select/deselect channels','Callback','spm_eeg_review_callbacks(''visu'',''sensor_select'',0)');
% handles.VIZU.sensors.select_all = uimenu(handles.VIZU.f_sensors,'Label','Select all channels','Callback','spm_eeg_review_callbacks(''visu'',''sensor_select_all'',0)');
% handles.VIZU.visu_misc = uimenu(f,'Label','Miscellaneous');
% handles.VIZU.misc.MainSwitch = uimenu(handles.VIZU.visu_misc,'Label','Reverse data sign','Callback','spm_eeg_review_callbacks(''visu'',''MainSwitch'',0)');
% handles.VIZU.misc.ygrid = uimenu(handles.VIZU.visu_misc,'Label','y-axis grid: on/off','Callback','spm_eeg_review_callbacks(''visu'',''ygrid'',0)');
% handles.VIZU.misc.xgrid = uimenu(handles.VIZU.visu_misc,'Label','x-axis grid: #seconds','Callback','spm_eeg_review_callbacks(''visu'',''xgrid'',0)');
% if D.Nsamples <= 200
%     set(handles.VIZU.time_w1,'enable','off')
% end

% % SELECTIONS uimenu
% g = uimenu('Label','Selections','enable',enab);
% handles.SELECT.select_add = uimenu(g,'Label','Adding selections');
% handles.SELECT.select_plus = uimenu(handles.SELECT.select_add,'Label','Add event to current selection',...
%     'callback','spm_eeg_review_callbacks(''select'',''add'')');
% handles.SELECT.select_all = uimenu(handles.SELECT.select_add,'Label','Select the entire loaded data',...
%     'callback','spm_eeg_review_callbacks(''select'',''all'')');
% handles.SELECT.select_remove = uimenu(g,'Label','Removing selections');
% handles.SELECT.select_minus = uimenu(handles.SELECT.select_remove,'Label',...
%     'Remove last selected event from current selection',...
%     'callback','spm_eeg_review_callbacks(''select'',''remove'')','enable','off');
% handles.SELECT.select_nothing = uimenu(handles.SELECT.select_remove,'Label','Clear all selections',...
%     'callback','spm_eeg_review_callbacks(''select'',''nothing'')','enable','off');
% handles.SELECT.select_view = uimenu(g,'Label','Viewing selections');
%
% handles.SELECT.goto_select1 = uimenu(handles.SELECT.select_view,'Label',...
%     'Go to closest selected event (forward)',...
%     'callback','spm_eeg_review_callbacks(''select'',''goto'',0)','enable','off');
% handles.SELECT.goto_select2 = uimenu(handles.SELECT.select_view,'Label',...
%     'Go to closest selected event (backward)',...
%     'callback','spm_eeg_review_callbacks(''select'',''goto'',1)','enable','off');
% handles.SELECT.show_select = uimenu(handles.SELECT.select_view,'Label',...
%     'Show concatenated events of current selection',...
%     'callback','spm_eeg_review_callbacks(''select'',''show'')','enable','off',...
%     'separator','on');
% handles.SELECT.Export_import_select = uimenu(g,'Label','Exporting/importing selections');
% handles.SELECT.save_select = uimenu(handles.SELECT.Export_import_select,'Label',...
%     'Save events of current selection',...
%     'callback','spm_eeg_review_callbacks(''select'',''save'')','enable','off');
% handles.SELECT.load_select = uimenu(handles.SELECT.Export_import_select,'Label',...
%     'Load selection from file',...
%     'callback','spm_eeg_review_callbacks(''select'',''load'')','enable','on');
%
% % TOOLS uimenu
% t = uimenu('Label','Tools','enable',enab);
% handles.TOOLS.event_av = uimenu(t,'Label','Event extraction');
% handles.TOOLS.cor_average = uimenu(handles.TOOLS.event_av,'Label','Coregister and average events',...
%     'callback','tools(''average_peak'')','enable','off');
% handles.TOOLS.find_peaks = uimenu(handles.TOOLS.event_av,'Label','Detect and average event',...
%     'callback','tools(''find_peaks'')','enable','off');
% handles.TOOLS.classify_peaks = uimenu(handles.TOOLS.event_av,'Label','Classify and average events',...
%     'callback','tools(''classify_peaks'')','enable','off');
% handles.TOOLS.freq_ana = uimenu(t,'Label','Frequency analysis');
% handles.TOOLS.spectrum_events = uimenu(handles.TOOLS.freq_ana,'Label','Frequency spectrum of iso-type events',...
%     'callback','tools(''spectrum_events'')','enable','off');
% handles.TOOLS.spectrum_comp = uimenu(handles.TOOLS.freq_ana,'Label','Concatenated events VS rest',...
%     'callback','tools(''spectrum_comp'')','enable','off');
% handles.TOOLS.IP = uimenu(t,'Label','Inverse Problem');
% handles.TOOLS.IP_template = uimenu(handles.TOOLS.IP,'Label','Template-based source reconstruction',...
%     'callback','tools(''IP_template'')','enable','off');
% handles.TOOLS.IP_indiv = uimenu(handles.TOOLS.IP,'Label','Forward modelling',...
%     'callback','tools(''IP_indiv'')','enable','off');

% BUTTONS
% xlim = D.PSD.VIZU.xlim;
% scrsz = get(0,'ScreenSize');


% load spm_eeg_review_buttons.mat


% % Visualization buttons
% handles.BUTTONS.vb1 = uicontrol(hfig,'Position',[0.14 0.92 0.05 0.04].*repmat(POS(3:4),1,2),'cdata',Y3,...
%     'Callback','spm_eeg_review_callbacks(''visu'',''iten_sc'',2)',...
%   'tooltipstring','Increase contrast (intensity rescaling)');
% set(handles.BUTTONS.vb1,'units','normalized');
%
% handles.BUTTONS.vb2 = uicontrol(hfig,'Position',[0.2 0.92 0.05 0.04].*repmat(POS(3:4),1,2),'cdata',Y4,...
%     'Callback','spm_eeg_review_callbacks(''visu'',''iten_sc'',0.5)',...
%   'tooltipstring','Decrease contrast (intensity rescaling)');
% set(handles.BUTTONS.vb2,'units','normalized');
%
%
% handles.BUTTONS.vb5 = uicontrol(hfig,'Position',[0.26 0.92 0.05 0.04].*repmat(POS(3:4),1,2),'cdata',Y7,...
%     'callback','spm_eeg_review_callbacks(''visu'',''zoom'',1)',...
%   'tooltipstring','Zoom in (mouse box)');
% set(handles.BUTTONS.vb5,'units','normalized');
%
% handles.BUTTONS.vb1 = uicontrol(hfig,'Position',[0.34 0.92 0.05 0.04].*repmat(POS(3:4),1,2),'cdata',Y8,...
%     'Callback','spm_eeg_review_callbacks(''visu'',''scalp_interp'',1)',...
%   'tooltipstring','scalp interpolation (image scalp data)');
% set(handles.BUTTONS.vb1,'units','normalized');


% switch D.PSD.type
%
%     case 'continuous'
%
%         % Selection buttons
%         Nevents = length(D.trials.events);
%         if Nevents >0
%             enab = 'on';
%         else
%             enab = 'off';
%         end
%
%         handles.BUTTONS.sb1 = uicontrol(hfig,'Position',[0.56 0.92 0.05 0.04].*repmat(POS(3:4),1,2),...
%             'cdata',Y9,'callback','spm_eeg_review_callbacks(''select'',''add'')',...
%             'tooltipstring','Add event to current selection (2 mouse clicks)');
%         set(handles.BUTTONS.sb1,'units','normalized');
%
%         handles.BUTTONS.sb2 = uicontrol(hfig,'Position',[0.42 0.92 0.05 0.04].*repmat(POS(3:4),1,2),...
%             'cdata',Y10,'callback','spm_eeg_review_callbacks(''select'',''goto'',0)',...
%           'tooltipstring','Go to closest selected event (forward)','enable',enab);
%         set(handles.BUTTONS.sb2,'units','normalized');
%
%         handles.BUTTONS.sb3 = uicontrol(hfig,'Position',[0.48 0.92 0.05 0.04].*repmat(POS(3:4),1,2),...
%             'cdata',Y11,'callback','spm_eeg_review_callbacks(''select'',''goto'',1)',...
%           'tooltipstring','Go to closest selected event (backward)','enable',enab);
%         set(handles.BUTTONS.sb3,'units','normalized');
%
%         clear Y1
%
%     case 'epoched'
%
%         % Event selection button
%         handles.BUTTONS.pop1 = uicontrol(hfig,'Position',[0.42 0.92 0.25 0.03].*repmat(POS(3:4),1,2),...
%             'style','popupmenu','string',D.PSD.trials.TrLabels,...
%             'callback','spm_eeg_review_callbacks(''select'',''switch'')');
%         set(handles.BUTTONS.pop1,'units','normalized')
%
% %         % Type of visualisation button
% %         handles.BUTTONS.pop2 = uicontrol(hfig,'Position',[0.82 0.96 0.15 0.03].*repmat(POS(3:4),1,2),...
% %             'style','popupmenu','string',{'standard','scalp'},...
% %             'callback','spm_eeg_review_callbacks(''visu'',''main'')');
% %         set(handles.BUTTONS.pop2,'units','normalized')
%
% end





% update handles
% handles.VIZU.menu_vizu = f;
% handles.SELECT.menu_select = g;
% handles.TOOLS.menu_tools = t;




%%
function [D] = PSD_initUD(D)
% function D = initUD(D)
% This function initializes the userdata structure required for the
% PreSelectdata toolbox.

D = struct(D);

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
        PSD.type = 'continuous';
        % after epoching
    case 'single'
        D.PSD.type = 'epoched';
        nTrials = length(D.trials);
        D.PSD.trials.TrLabels = cell(nTrials,1);
        for i = 1:nTrials
            D.PSD.trials.TrLabels{i} = ['Trial ',num2str(i),': ',D.trials(i).label];
        end
        D.PSD.trials.current = 1;
    case 'evoked'
        D.PSD.type = 'epoched';
        nTrials = length(D.trials);
        D.PSD.trials.TrLabels = cell(nTrials,1);
        for i = 1:nTrials
            D.PSD.trials.TrLabels{i} = ['Trial ',num2str(i),' (average of ',...
                num2str(D.trials(i).repl),' events): ',D.trials(i).label];
            D.trials(i).events = [];
        end
        D.PSD.trials.current = 1;
end

if strcmp(D.transform.ID,'time')

    % detect montage
    % ns = length(D.sensors);
    nc = length(D.channels);
    % Initialize visualization variables

    nt = D.Nsamples;
    Ieeg  = find(strcmp('EEG',{D.channels.type}));
    Imeg  = find(strcmp('MEG',{D.channels.type}));
    visuSensors = [Ieeg(:);Imeg(:)];
    M = sparse(length(visuSensors),nc);
    M(sub2ind(size(M),1:length(visuSensors),visuSensors(:)')) = 1;
    if isempty(visuSensors)
        warndlg('No EEG/MEG channels found!')
        visuSensors = [1:1:nc]';
    end
    decim = max([floor(size(D.data.y,2)./1e4),1]);
    data = D.data.y(visuSensors,1:decim:nt,:);
    scale = repmat(D.data.scale(visuSensors),[1 size(data,2) size(data,3)]);
    data = scale.*data;
    sd = std(data(:));
    % sd = mean(std(data,[],2));
    offset = [0:1:length(visuSensors)-1]'*sd/2;
    % [data,b,X] = detrendY(data);

    v_data                  = 0.25.*data +repmat(offset,[1 size(data,2) size(data,3)]);
    ma                      = max(v_data(:))+sd;
    mi                      = min(v_data(:))-sd;
    ylim                    = [mi ma];

    % D.PSD.tools.coreg          = 0;
    % D.PSD.tools.undo           = [];
    % D.PSD.tools.redo           = [];
    % D.PSD.tools.EventsAnalysis = [];

    D.PSD.VIZU.visu_scale      = 0.25;
    D.PSD.VIZU.FontSize        = 10;
    D.PSD.VIZU.visuSensors     = visuSensors;
    D.PSD.VIZU.visu_offset     = sd;
    D.PSD.VIZU.offset          = offset;
    D.PSD.VIZU.xlim            = [1,min([200,nt])];
    D.PSD.VIZU.ylim            = ylim;
    D.PSD.VIZU.ylim0           = ylim;
    % D.PSD.VIZU.detrend.b       = b;
    D.PSD.VIZU.figname         = 'main visualization window';
    D.PSD.VIZU.montage.M       = M;
    D.PSD.VIZU.montage.clab    = {D.channels(visuSensors).label};

else

    Ieeg  = find(strcmp('EEG',{D.channels.type}));
    Imeg  = find(strcmp('MEG',{D.channels.type}));
    visuSensors = [Ieeg(:);Imeg(:)];
    D.PSD.VIZU.visuSensors     = visuSensors;
    D.PSD.VIZU.montage.clab    = {D.channels(visuSensors).label};
    D.PSD.type = 'epoched';
    D.PSD.trials.current = 1;

end









