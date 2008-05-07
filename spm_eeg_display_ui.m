function Heeg = spm_eeg_display_ui(varargin)
% user interface for displaying EEG/MEG channel data.
% Heeg = spm_eeg_display_ui(varargin)
%
% optional argument:
% S         - struct
%    fields of S:
%     D       - MEEG object
%     Hfig    - Figure (or axes) to work in (Defaults to SPM graphics window)
%     rebuild - indicator variable: if defined spm_eeg_display_ui has been
%                                   called after channel selection
%
% output:
%     Heeg      - Handle of resulting figure
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_display_ui.m 1560 2008-05-07 12:18:58Z stefan $

if nargin == 1
    S = varargin{1};
end

if nargin == 0 || ~isfield(S, 'rebuild')
    try
        D = S.D;
    catch
        D = spm_select(1, '\.mat$', 'Select EEG mat file');
        try
            D = spm_eeg_load(D);
            sD = struct(D); % transform to struct for access to some fields
        catch
            error(sprintf('Trouble reading file %s', D));
        end
    end

    %     if D.ntrials == 1
    %         errordlg({'Continuous data cannot be displayed (yet).', 'Epoch first please.'});
    %         return;
    %     end
    %
    try
        % Use your own window
        F = S.Hfig;
    catch
        % use SPM graphics window
        F = findobj('Tag', 'Graphics');
        if isempty(F)
            F = spm_figure('create','Graphics','Graphics','on');
        end
    end

    set(F, 'SelectionType', 'normal');

    handles = guihandles(F);

    handles.colour = {[0 0 1], [1 0 0], [0 1 0], [1 0 1]};
    handles.badchannels = D.badchannels;

    % variable needed to store current trial listbox selection
    handles.Tselection = 1;

    % fontsize used troughout
    % (better compute that fontsize)
    FS1 = spm('FontSize', 8);

    figure(F);clf

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % setup of GUI elements
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Slider for trial number
    %------------------------
    if D.ntrials > 1

        % text above slider
        uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
            'String', 'Trial number',...
            'Position',[0.11 0.072 0.15 0.029],...
            'HorizontalAlignment', 'right', 'FontSize', FS1,...
            'BackgroundColor', 'w');

        handles.trialslider = uicontrol(F, 'Tag', 'trialslider', 'Style', 'slider',...
            'Min', 1, 'Max', D.ntrials, 'Value', 1, 'Units',...
            'normalized', 'Position', [0.05 0.045 0.25 0.03],...
            'SliderStep', [1/(D.ntrials-1) min(D.ntrials-1, 10/(D.ntrials-1))],...
            'TooltipString', 'Choose trial number',...
            'Callback', @trialslider_update,...
            'Parent', F, 'Interruptible', 'off');

        % frame for trialslider text
        uicontrol(F, 'Style','Frame','BackgroundColor', spm('Colour'), 'Units',...
            'normalized', 'Position',[0.05 0.019 0.25 0.031]);

        % trials slider texts
        uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
            'String', '1',...
            'Position',[0.06 0.02 0.07 0.029],...
            'HorizontalAlignment', 'left', 'FontSize', FS1,...
            'BackgroundColor', spm('Colour'));

        handles.trialtext = uicontrol(F, 'Style', 'text', 'Tag', 'trialtext',...
            'Units', 'normalized',...
            'String', int2str(get(handles.trialslider, 'Value')),...
            'Position',[0.14 0.02 0.07 0.029],...
            'HorizontalAlignment', 'center', 'FontSize', FS1,...
            'BackgroundColor', spm('Colour'));

        uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
            'String', mat2str(D.ntrials),...
            'Position',[0.23 0.02 0.06 0.029],...
            'HorizontalAlignment', 'right', 'FontSize', FS1,...
            'BackgroundColor', spm('Colour'));

    end

    % Scaling slider
    %-----------------
    % estimate of maximum scaling value
    if strcmp(D.transformtype, 'TF')
        handles.scalemax = 2*max(max(max(abs(D(setdiff(D.meegchannels, D.badchannels), :,:, 1)))));
    else
        handles.scalemax = 2*max(max(abs(D(setdiff(D.meegchannels, D.badchannels), :, 1))));
    end

    handles.order = 10^floor(log10(handles.scalemax)-2);
    scale = 100;

    % text above slider
    uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
        'String', 'Scaling',...
        'Position',[0.37 0.072 0.15 0.029],...
        'HorizontalAlignment', 'right', 'FontSize', FS1,...
        'BackgroundColor', 'w');

    % slider
    handles.scaleslider = uicontrol('Tag', 'scaleslider', 'Style', 'slider',...
        'Min', 1, 'Max', 300, 'Value', scale, 'Units',...
        'normalized', 'Position', [0.35 0.045 0.25 0.03],...
        'SliderStep', [0.0033 0.1],...
        'TooltipString', 'Choose scaling',...
        'Callback', @scaleslider_update,...
        'Parent', F);

    % frame for text
    uicontrol(F, 'Style','Frame','BackgroundColor', spm('Colour'), 'Units',...
        'normalized', 'Position',[0.35 0.019 0.25 0.031]);

    % text
    uicontrol(F, 'Style', 'text', 'Units', 'normalized', ...
        'String', num2str(handles.order, 4),...
        'Position',[0.36 0.02 0.07 0.029],...
        'HorizontalAlignment', 'left', 'FontSize', FS1,...
        'BackgroundColor', spm('Colour'));

    handles.scaletext = uicontrol(F, 'Style', 'text', 'Tag', 'scaletext',...
        'Units', 'normalized',...
        'String', num2str(round(get(handles.scaleslider, 'Value'))*handles.order, 4),...
        'Position',[0.44 0.02 0.07 0.029],...
        'HorizontalAlignment', 'center', 'FontSize', FS1,...
        'BackgroundColor',spm('Colour'));

    uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
        'String', num2str(round(get(handles.scaleslider, 'Max'))*handles.order, 4),...
        'Position',[0.52 0.02 0.07 0.029],...
        'HorizontalAlignment', 'right', 'FontSize', FS1,...
        'BackgroundColor',spm('Colour'));


    %---------------------
    % Save pushbutton
    uicontrol('Tag', 'savebutton', 'Style', 'pushbutton',...
        'Units', 'normalized', 'Position', [0.615 0.02 0.11 0.03],...
        'String', 'Save', 'FontSize', FS1,...
        'BackgroundColor', spm('Colour'),...
        'CallBack', @savebutton_update,...
        'Parent', F);

    % Reject selection pushbutton
    uicontrol('Tag', 'channelselectbutton', 'Style', 'pushbutton',...
        'Units', 'normalized', 'Position', [0.615 0.05 0.11 0.03],...
        'String', 'Reject', 'FontSize', FS1,...
        'BackgroundColor', spm('Colour'),...
        'CallBack', @rejectbutton_update,...
        'Parent', F);

    % Channel selection pushbutton
    uicontrol('Tag', 'channelselectbutton', 'Style', 'pushbutton',...
        'Units', 'normalized', 'Position', [0.615 0.08 0.11 0.03],...
        'String', 'Channels', 'FontSize', FS1,...
        'BackgroundColor', spm('Colour'),...
        'CallBack', @channelselectbutton_update,...
        'Parent', F);

    % Topography display pushbutton
    uicontrol('Tag', '3Dtopographybutton', 'Style', 'pushbutton',...
        'Units', 'normalized', 'Position', [0.615 0.11 0.11 0.03],...
        'String', 'Topography', 'FontSize', FS1,...
        'BackgroundColor', spm('Colour'),...
        'CallBack', @scalp3d_select,...
        'Parent', F);

    % trial listbox
    if D.ntrials > 1
        trialnames = {};
        for i = 1:D.ntrials
            trialnames{i} = [sprintf('%-12s', sprintf('trial %d', i))  sprintf('%-4s', strvcat(D.conditions(i)))];
            if D.reject(i)
                trialnames{i} = [trialnames{i} ' reject'];
            end
        end
        handles.trialnames = trialnames;

        handles.triallistbox = uicontrol('Tag', 'triallistbox', 'Style', 'listbox',...
            'Units', 'normalized', 'Position', [0.74 0.02 0.2 0.2],...
            'Min', 0, 'Max', 2, 'String', trialnames,...
            'HorizontalAlignment', 'left',...
            'Value', 1,...
            'BackgroundColor', [1 1 1],...
            'CallBack', @triallistbox_update,...
            'Parent', F);
    end

    % display file name at top of page
    uicontrol('Style', 'text', 'String', fullfile(D.path, D.fname),...
        'Units', 'normalized', 'Position', [0.05 0.98 0.8 0.02],...
        'Background', 'white', 'FontSize', 14, 'HorizontalAlignment', 'Left');

    % axes with scaling and ms display
    axes('Position', [0.05 0.15 0.2 0.07], 'YLim', [-round(scale)*handles.order round(scale)*handles.order], 'XLim', [1 D.nsamples],...
        'XTick', [], 'YTick', [], 'LineWidth', 2);
    handles.scaletexttop = text('Position', [0, round(scale)*handles.order], 'String', ['   ' num2str(round(scale)*handles.order, 4)], 'Interpreter', 'Tex',...
        'FontSize', FS1, 'VerticalAlignment', 'top',...
        'Tag', 'scaletext2');
    try
        ylabel(D.units(1), 'Interpreter', 'tex', 'FontSize', FS1);
    end
    text('Position', [D.nsamples, -round(scale)*handles.order], 'String', sprintf('%d', round((D.nsamples-1)*1000/D.fsample)), 'Interpreter', 'Tex',...
        'FontSize', FS1, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    xlabel('ms', 'Interpreter', 'tex', 'FontSize', FS1);

    % Added option to specify in advance only subset of channels to display RH
    % ([] = prompt for channel file; char = load channel file; No error checking yet)
    try
        S.gfx.channels = S.chans;
        if isempty(S.chans)
            S.chans = spm_select(1, '\.mat$', 'Select channel mat file');
            S.load = load(S.chans);
            S.gfx.channels = S.load.Iselectedchannels;
        elseif ischar(S.chans)
            S.load = load(S.chans);
            S.gfx.channels = S.load.Iselectedchannels;
        end
    catch
        % channels to display, initially exclude bad channels
        S.gfx.channels = setdiff(D.meegchannels, D.badchannels);
    end
    gfx = S.gfx; D = cache(D,gfx);

else
    % this is a re-display with different set of selected channels, delete plots
    handles = guidata(S.Hfig);
    delete(handles.Heegaxes);
    handles.Heegaxes = [];

    for i = 1:length(handles.Heegfigures)
        if ~isempty(handles.Heegfigures{i})
            delete(handles.Heegfigures{i});
        end
    end

    D = S.D;
    F = S.Hfig;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplots of EEG data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% position of plotting area for eeg data in graphics figure
Pos = [0.03 0.23 0.95 0.75];

% Compute width of display boxes
%-------------------------------

% indices of displayed channels (in order of data)
handles.Cselection = S.gfx.channels;
p = coor2D(D, handles.Cselection);
Rxy = 1.5; % ratio of x- to y-axis lengths

Npos = size(p, 2); % number of positions

if Npos > 1
    % more than 1 channel for display
    for i = 1:Npos
        for j = 1:Npos
            % distance between channels
            d(i,j) = sqrt(sum((p(:,j)-p(:,i)).^2));

            % their angle
            alpha(i,j) = acos((p(1,j)-p(1,i))/(d(i,j)+eps));
        end
    end
    d = d/2;

    alpha(alpha > pi/2) = pi-alpha(alpha > pi/2);
    Talpha = asin(1/(sqrt(1+Rxy^2)));

    for i = 1:Npos
        for j = 1:Npos
            if alpha(i,j) <= Talpha
                x(i,j) = d(i,j)*cos(alpha(i,j));
            else
                x(i,j) = Rxy*d(i,j)*cos(pi/2-alpha(i,j));
            end
        end
    end


    % half length of axes in x-direction
    Lxrec = min(x(find(x~=0)));

else
    % only one channel
    Lxrec = 1;
end

% coordinates of lower left corner of drawing boxes
p(1, :) = p(1, :) - Lxrec;
p(2, :) = p(2, :) - Lxrec/Rxy;

% envelope of coordinates
e = [min(p(1,:)) max(p(1,:))+2*Lxrec min(p(2,:)) max(p(2,:))+2*Lxrec/Rxy];

% shift coordinates to zero
p(1,:) = p(1,:) - mean(e(1:2));
p(2,:) = p(2,:) - mean(e(3:4));

% scale such that envelope goes from -0.5 to 0.5
Sf = 0.5/max(max(abs(p(1,:))), (max(abs(p(2,:)))));
p = Sf*p;
Lxrec = Sf*Lxrec;

% and back to centre
p = p+0.5;

% translate and scale to fit into drawing area of figure
p(1,:) = Pos(3)*p(1,:)+Pos(1);
p(2,:) = Pos(4)*p(2,:)+Pos(2);

scale = get(handles.scaleslider, 'Value')*handles.order;

% cell vector for figures handles of separate single channel plots
handles.Heegfigures = cell(1, Npos);

% cell vector for axes handles of single channel plots
handles.Heegaxes2 = cell(1, Npos);

% plot the graphs
for i = 1:Npos

    % uicontextmenus for axes
    Heegmenus(i) = uicontextmenu;
    if ismember(i, D.badchannels)
        labelstring = 'Declare as good';
    else
        labelstring = 'Declare as bad';
    end

    uimenu(Heegmenus(i), 'Label',...
        sprintf('%s (%d)', strvcat(D.chanlabels(handles.Cselection(i))), handles.Cselection(i)));
    uimenu(Heegmenus(i), 'Separator', 'on');
    uimenu(Heegmenus(i), 'Label', labelstring,...
        'CallBack', {@switch_bad, i});

    handles.Heegaxes(i) = axes('Position',...
        [p(1,i) p(2,i) 2*Lxrec*Pos(3) 2*Lxrec/Rxy*Pos(4)],...
        'ButtonDownFcn', {@windowplot, i},...
        'NextPlot', 'add',...
        'Parent', F,...
        'UIContextMenu', Heegmenus(i),...
        'YLim', [-scale scale],...
        'XLim', [D.time(1) D.time(end)],...
        'XTick', [], 'YTick', [], 'Box', 'off');

    for j = 1:length(handles.Tselection)

        if strcmp(D.transformtype, 'TF')
            h = imagesc(squeeze(D(handles.Cselection(i), :,:, handles.Tselection(j))));
            set(h, 'ButtonDownFcn', {@windowplot, i},...
                'Clipping', 'off', 'UIContextMenu', Heegmenus(i));
        else
            h = plot(handles.Heegaxes(i), D.time,...
                D(handles.Cselection(i), :, handles.Tselection(j)),...
                'Color', handles.colour{j},...
                'ButtonDownFcn', {@windowplot, i},...
                'Clipping', 'off', 'UIContextMenu', Heegmenus(i));
        end

        if strcmp(D.transformtype, 'TF')
            set(gca, 'ZLim', [-scale scale],...
                'XLim', [1 D.nsamples], 'YLim',  [1 D.nfrequencies], 'XTick', [], 'YTick', [], 'ZTick', [],'Box', 'off');
            caxis([-scale scale])
            colormap('jet')
        else
            if Lxrec > 0.1
                % boxes are quite large
                set(handles.Heegaxes(i), 'Xgrid', 'on');
            end
        end
    end
end

handles.Lxrec = Lxrec;
handles.D = D;

% store handles
guidata(F, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% callbacks for GUI elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function switch_bad(hObject, events, ind)
% update called from channel contextmenu

handles = guidata(hObject);

ind = handles.Cselection(ind);

if ismember(ind, handles.badchannels)
    handles.badchannels = setdiff(handles.badchannels, ind);
else
    handles.badchannels = unique([handles.badchannels, ind]);
end

Heegmenu = uicontextmenu;
if ismember(ind, handles.badchannels)
    labelstring = sprintf('Declare %s as good', handles.D.chanlabels(ind));
else
    labelstring = sprintf('Declare %s as bad', handles.D.chanlabels(ind));
end

uimenu(Heegmenu, 'Label', labelstring,...
    'CallBack', {@switch_bad, ind});

set(handles.Heegaxes(ind), 'UIContextMenu', Heegmenu);

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trialslider_update(hObject, events)
% update called from trialslider

handles = guidata(hObject);

% slider value
ind = round(get(hObject, 'Value'));
set(hObject, 'Value', ind);

% update triallistbox
set(handles.triallistbox, 'Value', ind);

% Update display of current trial number
set(handles.trialtext, 'String', mat2str(ind));

% make plots
draw_subplots(handles, ind)

handles.Tselection = ind;
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function triallistbox_update(hObject, events)
% update called from triallistbox

handles = guidata(hObject);

% listbox selection
ind = get(hObject, 'Value');

if length(ind) > length(handles.colour)
    warning('Can only display %d different traces', length(handles.colour));
    ind = handles.Tselection;
    set(hObject, 'Value', handles.Tselection);
elseif length(ind) < 1
    % can happen if user presses with cntl on already selected trial
    ind = handles.Tselection;
    set(hObject, 'Value', handles.Tselection);
else
    % update trialslider to minimum of selection
    set(handles.trialslider, 'Value', min(ind));

    % Update display of current trial number
    set(handles.trialtext, 'String', mat2str(min(ind)));

    handles.Tselection = ind;
    guidata(hObject, handles);

    % make plots
    draw_subplots(handles, ind)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scaleslider_update(hObject, events)
% update called from scaleslider

handles = guidata(hObject);
D = handles.D;

% slider value
scale = round(get(hObject, 'Value'))*handles.order;
% set(hObject, 'Value', scale);
%
% text below slider
set(handles.scaletext, 'String', num2str(scale, 4));

% text at top of figure
set(handles.scaletexttop, 'String', ['   ' num2str(scale, 4)]);

% rescale plots
for i = 1:length(handles.Heegaxes)

    
    if strcmp(D.transformtype, 'TF')
        set(handles.Heegaxes(i), 'ZLim', [0 scale],...
            'XLim', [1 D.nsamples], 'YLim', [1 D.nfrequencies], 'XTick', [], 'YTick', [], 'ZTick', [],'Box', 'off');
        caxis([-scale scale])
    else
        if handles.Lxrec > 0.1
            % boxes are quite large
            set(handles.Heegaxes(i), 'YLim', [-scale scale],...
                'XLim',[D.time(1) D.time(end)], 'Box', 'off', 'Xgrid', 'on');
        else
            % otherwise remove tickmarks
            set(handles.Heegaxes(i), 'YLim', [-scale scale],...
                'XLim', [D.time(1) D.time(end)], 'XTick', [], 'YTick', [], 'Box', 'off');
        end
    end
end

% rescale separate windows (if there are any)
for i = 1:length(handles.Heegfigures)

    if ~isempty(handles.Heegfigures{i})

        if strcmp(D.transformtype, 'TF')
            caxis(handles.Heegaxes2{i}, [-scale scale])
        else
            set(handles.Heegaxes2{i}, 'YLim', [-scale scale],...
                'XLim', [1000*D.time(1) 1000*D.time(end)]);
        end

        % update legend
        legend;

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function draw_subplots(handles, ind)
% This function plots data in the subplots

D = handles.D;

scale = round(get(handles.scaleslider, 'Value'));

for i = 1:length(handles.Heegaxes)

    cla(handles.Heegaxes(i));

    for j = 1:length(ind)
        if strcmp(D.transformtype, 'TF')

            h = imagesc(squeeze(D(handles.Cselection(i), :,:, ind(j))), 'Parent', handles.Heegaxes(i));
            set(h, 'ButtonDownFcn', {@windowplot, i},...
                'Clipping', 'off');
        else
            h = plot(handles.Heegaxes(i), D.time,...
                D(handles.Cselection(i), :, ind(j)),...
                'Color', handles.colour{j},'ButtonDownFcn', {@windowplot, i},...
                'Clipping', 'off');
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function windowplot(hObject, events, ind)
% this function plots data from channel ind in a separate window

st = get(gcf,'SelectionType');

% do nothing for right button click
if strcmp(st, 'alt')
    return;
end

handles = guidata(hObject);
D = handles.D;

w = handles.Heegfigures{ind};

% index to get the channel name
if ~isempty(w)

    % delete the existing window
    delete(w);
    handles.Heegaxes2{ind} = [];
    handles.Heegfigures{ind} = [];
else
    handles.Heegfigures{ind} = figure;
    set(gcf,...
        'Name', (sprintf('Channel %s', strvcat(D.chanlabels(handles.Cselection(ind))))),...
        'NumberTitle', 'off',...
        'DeleteFcn', {@delete_Heegwindows, ind});

    handles.Heegaxes2{ind} = gca;
    set(handles.Heegaxes2{ind}, 'NextPlot', 'add');

    if strcmp(D.transformtype, 'TF')
        xlabel('ms', 'FontSize', 16);
        ylabel('Hz', 'FontSize', 16, 'Interpreter', 'Tex')

        for i = 1:length(handles.Tselection)

            imagesc(1000*D.time, D.frequencies, squeeze(D(handles.Cselection(ind), :,:, handles.Tselection(i))));

        end

        scale = get(handles.scaleslider, 'Value');
        set(handles.Heegaxes2{ind}, 'ZLim', [-scale scale],...
            'XLim',  1000*[D.time(1) D.time(end)], 'YLim', [min(D.frequencies) max(D.frequencies)], 'Box', 'on');
        colormap('jet')
        caxis([-scale scale])

    else
        xlabel('ms', 'FontSize', 16);
        ylabel(D.units(handles.Cselection(ind)), 'FontSize', 16, 'Interpreter', 'Tex')

        for i = 1:length(handles.Tselection)
            plot(1000*D.time,...
                D(handles.Cselection(ind), :, handles.Tselection(i)),...
                'Color', handles.colour{i}, 'LineWidth', 2);
        end
        set(gca, 'YLim', [-scale scale],...
            'XLim', 1000*[D.time(1) D.time(end)], 'Box', 'on');
        grid on
        
        if D.ntrials > 1
            legend(handles.trialnames{handles.Tselection}, 0);
        end

    
    end

    scale = get(handles.scaleslider, 'Value')*handles.order;


    title(sprintf('%s (%d)', strvcat(D.chanlabels(handles.Cselection(ind))), handles.Cselection(ind)), 'FontSize', 16);


    % save handles structure to new figure
    guidata(handles.Heegaxes2{ind}, handles);
end

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function delete_Heegwindows(hObject, events, ind)
% deletes window of ind-th channel when delete was not via the eeg channel
% plot in the SPM graphics window

% returns handles structure of the graphics window
handles = guidata(hObject);
handles = guidata(handles.Graphics);

delete(handles.Heegfigures{ind});
handles.Heegfigures{ind} = [];

% save to handles structure in main Graphics figure
guidata(handles.Graphics, handles);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function savebutton_update(hObject, events)
% save updated matfile
handles = guidata(hObject);
D = handles.D;

spm('Pointer', 'Watch');

save(D);

spm('Pointer', 'Arrow');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rejectbutton_update(hObject, events)
% called from reject button

handles = guidata(hObject);
D = handles.D;

% listbox selection
ind = handles.Tselection;

if ind > 1
    % toggle first select only
    ind = ind(1);
end

if D.reject(ind) == 0
    % reject this trial
    tmp = [sprintf('%-12s', sprintf('trial %d', ind))  sprintf('%s', D.conditions(ind))...
        sprintf(' %s', 'reject')];
    D = reject(D, ind, 1);
else
    % un-reject this trial
    tmp = [sprintf('%-12s', sprintf('trial %d', ind))  sprintf('%s', D.conditions(ind))];
    D = reject(D, ind, 0);
end
handles.trialnames{ind} = tmp;

set(handles.triallistbox, 'String', handles.trialnames);
handles.D = D;

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function channelselectbutton_update(hObject, events)
handles = guidata(hObject);
D = handles.D;

ind = spm_eeg_select_channels(D);
gfx.channels = ind;

D = cache(D, gfx);
S.D = D;
S.gfx = gfx;
S.rebuild = 1;
S.Hfig = handles.Graphics;

spm_eeg_display_ui(S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scalp3d_select(hObject, events)
% ask user for peri-stimulus time and call scalp3d

handles = guidata(hObject);
D = handles.D;

ok = 0;
while ~ok
    try
        answer = spm_eeg_scalp_dlg;
    catch
        return
    end
    t = answer{1};
    s = answer{2};
    for n=1:length(t)
        if t(n) >= D.time(1,'ms')...
                && t(n) <= D.time(D.nsamples, 'ms')...
                && (strcmpi(s, '2d') || strcmpi(s, '3d'))
            ok = 1;
        end
    end
end

% call scalp2d or 3d
%--------------------------------------------------------------------------
spm('Pointer', 'Watch');drawnow;
if strcmpi(s, '2d')
    spm_eeg_scalp2d_ext(D, t, handles.Tselection(1));
else
    %     % change time for james' function
    %     t=round(t/1000*D.Radc)+D.events.start+1;
    %     if length(t) == 1
    %         d = squeeze(D.data(D.channels.eeg, t, handles.Tselection(1)));
    %     else
    %         d = squeeze(mean(D.data(D.channels.eeg, t, handles.Tselection(1)), 2));
    %     end
    %     if strmatch(D.channels.ctf,'bdf_setup.mat')
    %         spm_eeg_scalp3d(d);
    %     else
    %         errordlg({'This can only be used for 128 channel BDF data at the moment'});
    %         return;
    %     end
    %

end
spm('Pointer', 'Arrow');drawnow;

