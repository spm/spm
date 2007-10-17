function Heeg = spm_eeg_display_ui(varargin)
% user interface for displaying EEG/MEG channel data.
% Heeg = spm_eeg_display_ui(varargin)
%
% optional argument:
% S		    - struct
%    fields of S:
%     D 	  - EEG struct
%     Hfig	  - Figure (or axes) to work in (Defaults to SPM graphics window)
%     rebuild - indicator variable: if defined spm_eeg_display_ui has been
%                                   called after channel selection
%
% output:
%     Heeg		- Handle of resulting figure
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_display_ui.m 955 2007-10-17 15:15:09Z rik $

if nargin == 1
    S = varargin{1};
end

if nargin == 0 | ~isfield(S, 'rebuild')
    try
        D = S.D;
    catch
        D = spm_select(1, '\.mat$', 'Select EEG mat file');
        try
            D = spm_eeg_ldata(D);
        catch
            error(sprintf('Trouble reading file %s', D));
        end
    end

    if ~isfield(D.channels, 'Bad')
        D.channels.Bad = [];
    end


    if D.Nevents == 1 & ~isfield(D.events, 'start')
        errordlg({'Continuous data cannot be displayed (yet).', 'Epoch first please.'});
        return;
    end

    % units, default EEG
    if ~isfield(D, 'units')
        D.units = '\muV';
    end

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

    % variable needed to store current trial listbox selection
    handles.Tselection = 1;

    % fontsize used troughout
    % (better compute that fontsize)
%    FS1 = spm('FontSize', 14);
    FS1 = spm('FontSize', 8);

    figure(F);clf

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % setup of GUI elements
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Slider for trial number
    %------------------------
    if D.Nevents > 1

        % text above slider
        uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
            'String', 'Trial number',...
            'Position',[0.11 0.072 0.15 0.029],...
            'HorizontalAlignment', 'right', 'FontSize', FS1,...
            'BackgroundColor', 'w');

        handles.trialslider = uicontrol(F, 'Tag', 'trialslider', 'Style', 'slider',...
            'Min', 1, 'Max', D.Nevents, 'Value', 1, 'Units',...
            'normalized', 'Position', [0.05 0.045 0.25 0.03],...
            'SliderStep', [1/(D.Nevents-1) min(D.Nevents-1, 10/(D.Nevents-1))],...
            'TooltipString', 'Choose trial number',...
            'Callback', @trialslider_update,...
            'Parent', F, 'Interruptible', 'off');

        % frame for trialslider text
        uicontrol(F, 'Style','Frame','BackgroundColor',spm('Colour'), 'Units',...
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
            'String', mat2str(D.Nevents),...
            'Position',[0.23 0.02 0.06 0.029],...
            'HorizontalAlignment', 'right', 'FontSize', FS1,...
            'BackgroundColor', spm('Colour'));

    end

    % Scaling slider
    %-----------------
    % estimate of maximum scaling value
    if isfield (D,'Nfrequencies')
        handles.scalemax = 2*ceil(max(max(max(abs(D.data(setdiff(D.channels.eeg, D.channels.Bad), :,:, 1))))));
    else
%        handles.scalemax = 2*ceil(max(max(max(abs(D.data(setdiff([1:D.Nchannels], D.channels.Bad), :, :))))));
        handles.scalemax = 2*ceil(max(max(abs(D.data(setdiff(D.channels.eeg, D.channels.Bad), :, 1)))));
    end
    scale = handles.scalemax/2;

    % text above slider
    uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
        'String', 'Scaling',...
        'Position',[0.37 0.072 0.15 0.029],...
        'HorizontalAlignment', 'right', 'FontSize', FS1,...
        'BackgroundColor', 'w');

    % slider
    handles.scaleslider = uicontrol('Tag', 'scaleslider', 'Style', 'slider',...
        'Min', 1, 'Max', handles.scalemax, 'Value', scale, 'Units',...
        'normalized', 'Position', [0.35 0.045 0.25 0.03],...
        'SliderStep', [1/(handles.scalemax-1) 10/(handles.scalemax-1)],...
        'TooltipString', 'Choose scaling',...
        'Callback', @scaleslider_update,...
        'Parent', F);

    % frame for text
    uicontrol(F, 'Style','Frame','BackgroundColor',spm('Colour'), 'Units',...
        'normalized', 'Position',[0.35 0.019 0.25 0.031]);

    % text
    uicontrol(F, 'Style', 'text', 'Units', 'normalized', 'String', '1',...
        'Position',[0.36 0.02 0.07 0.029],...
        'HorizontalAlignment', 'left', 'FontSize', FS1,...
        'BackgroundColor',spm('Colour'));

    handles.scaletext = uicontrol(F, 'Style', 'text', 'Tag', 'scaletext',...
        'Units', 'normalized',...
        'String', mat2str(handles.scalemax/2),...
        'Position',[0.44 0.02 0.07 0.029],...
        'HorizontalAlignment', 'center', 'FontSize', FS1,...
        'BackgroundColor',spm('Colour'));

    uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
        'String', mat2str(handles.scalemax),...
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
    if D.Nevents > 1
        trialnames = {};
        for i = 1:D.Nevents
            trialnames{i} = [sprintf('%-12s', sprintf('trial %d', i))  sprintf('%-4s', sprintf('%d', D.events.code(i)))];
            if D.events.reject(i)
                trialnames{i} = [trialnames{i} sprintf('%-8s', 'reject')];
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
    axes('Position', [0.05 0.15 0.2 0.07]);
    set(gca, 'YLim', [-scale scale], 'XLim', [1 D.Nsamples],...
        'XTick', [], 'YTick', [], 'LineWidth', 2);
    handles.scaletexttop = text(0, scale, sprintf(' %d', 2*scale), 'Interpreter', 'Tex',...
        'FontSize', FS1, 'VerticalAlignment', 'top',...
        'Tag', 'scaletext2');
    ylabel(D.units, 'Interpreter', 'tex', 'FontSize', FS1);
    text(D.Nsamples, -scale, sprintf('%d', round((D.Nsamples-1)*1000/D.Radc)), 'Interpreter', 'Tex',...
        'FontSize', FS1, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    xlabel('ms', 'Interpreter', 'tex', 'FontSize', FS1);

    % Added option to specify in advance only subset of channels to display RH
    % ([] = prompt for channel file; char = load channel file; No error checking yet)
    try
    	D.gfx.channels = S.chans;
        if isempty(S.chans)
            S.chans = spm_select(1, '\.mat$', 'Select channel mat file');
	    S.load = load(S.chans);
    	    D.gfx.channels = S.load.Iselectedchannels;
        elseif ischar(S.chans)
	    S.load = load(S.chans);
    	    D.gfx.channels = S.load.Iselectedchannels;
	end
    catch
    % channels to display, initially exclude bad channels
        D.gfx.channels = setdiff([1:length(D.channels.order)], D.channels.Bad);
    end

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
Csetup = load(fullfile(spm('dir'), 'EEGtemplates', D.channels.ctf));

% indices of displayed channels (in order of data)
handles.Cselection2 = D.gfx.channels;

% indices of displayed channels (in order of the channel template file)
handles.Cselection = D.channels.order(handles.Cselection2);

p = Csetup.Cpos(:, handles.Cselection);
Rxy = Csetup.Rxy; % ratio of x- to y-axis lengths

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

scale = get(handles.scaleslider, 'Value');

% cell vector for figures handles of separate single channel plots
handles.Heegfigures = cell(1, Npos);

% cell vector for axes handles of single channel plots
handles.Heegaxes2 = cell(1, Npos);

% plot the graphs
for i = 1:Npos

    % uicontextmenus for axes
    Heegmenus(i) = uicontextmenu;
    if ismember(i, D.channels.Bad)
        labelstring = 'Declare as good';
    else
        labelstring = 'Declare as bad';
    end

    uimenu(Heegmenus(i), 'Label',...
        sprintf('%s (%d)', D.channels.name{handles.Cselection2(i)}, handles.Cselection2(i)));
    uimenu(Heegmenus(i), 'Separator', 'on');
    uimenu(Heegmenus(i), 'Label', labelstring,...
        'CallBack', {@switch_bad, i});

    handles.Heegaxes(i) = axes('Position',...
        [p(1,i) p(2,i) 2*Lxrec*Pos(3) 2*Lxrec/Rxy*Pos(4)],...
        'ButtonDownFcn', {@windowplot, i},...
        'NextPlot', 'add',...
        'Parent', F,...
        'UIContextMenu', Heegmenus(i));

    % make axes current
    axes(handles.Heegaxes(i));

    for j = 1:length(handles.Tselection)

        if isfield(D,'Nfrequencies')
            h = imagesc(squeeze(D.data(handles.Cselection2(i), :,:, handles.Tselection(j))));
            set(h, 'ButtonDownFcn', {@windowplot, i},...
                'Clipping', 'off', 'UIContextMenu', Heegmenus(i));
        else
            h = plot([-D.events.start:D.events.stop]*1000/D.Radc,...
                D.data(handles.Cselection2(i), :, handles.Tselection(j)),...
                'Color', handles.colour{j});
            set(h, 'ButtonDownFcn', {@windowplot, i},...
                'Clipping', 'off', 'UIContextMenu', Heegmenus(i));
        end
    end
    if isfield(D,'Nfrequencies')
        set(gca, 'ZLim', [-scale scale],...
            'XLim', [1 D.Nsamples], 'YLim',  [1 D.Nfrequencies], 'XTick', [], 'YTick', [], 'ZTick', [],'Box', 'off');
        caxis([-scale scale])
        colormap('jet')
    else
        if Lxrec > 0.1
            % boxes are quite large
            set(gca, 'YLim', [-scale scale],...
                'XLim', [-D.events.start D.events.stop]*1000/D.Radc, 'Box', 'off', 'Xgrid', 'on');
        else
            % otherwise remove tickmarks
            set(gca, 'YLim', [-scale scale],...
                'XLim', [-D.events.start D.events.stop]*1000/D.Radc, 'XTick', [], 'YTick', [], 'Box', 'off');
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

% ind1 = handles.Cselection(ind);
ind2 = handles.Cselection2(ind);

if ismember(ind2, handles.D.channels.Bad)
    handles.D.channels.Bad = setdiff(handles.D.channels.Bad, ind2);
else
    handles.D.channels.Bad = unique([handles.D.channels.Bad ind2]);
end

Heegmenu = uicontextmenu;
if ismember(ind, handles.D.channels.Bad)
    labelstring = sprintf('Declare %s as good', handles.D.channels.name{ind2});
else
    labelstring = sprintf('Declare %s as bad', handles.D.channels.name{ind2});
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
scale = round(get(hObject, 'Value'));
set(hObject, 'Value', scale);

% text below slider
set(handles.scaletext, 'String', mat2str(scale));

% text at top of figure
set(handles.scaletexttop, 'String', sprintf(' %d', 2*scale));

% rescale plots
for i = 1:length(handles.Heegaxes)

    % make ith subplot current
    axes(handles.Heegaxes(i));

    if isfield(D,'Nfrequencies')
        set(gca, 'ZLim', [0 scale],...
            'XLim', [1 D.Nsamples], 'YLim', [1 D.Nfrequencies], 'XTick', [], 'YTick', [], 'ZTick', [],'Box', 'off');
        caxis([-scale scale])
    else
        if handles.Lxrec > 0.1
            % boxes are quite large
            set(gca, 'YLim', [-scale scale],...
                'XLim', [-D.events.start D.events.stop]*1000/D.Radc, 'Box', 'off', 'Xgrid', 'on');
        else
            % otherwise remove tickmarks
            set(gca, 'YLim', [-scale scale],...
                'XLim', [-D.events.start D.events.stop]*1000/D.Radc, 'XTick', [], 'YTick', [], 'Box', 'off');
        end
    end
end

% rescale separate windows (if there are any)
for i = 1:length(handles.Heegfigures)
    if ~isempty(handles.Heegfigures{i})

        axes(handles.Heegaxes2{i});

        if isfield(D,'Nfrequencies')

            caxis([-scale scale])
        else
            set(gca, 'YLim', [-scale scale],...
                'XLim', [-D.events.start D.events.stop]*1000/D.Radc);
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

    % make ith subplot current
    axes(handles.Heegaxes(i));

    cla
    set(handles.Heegaxes(i), 'NextPlot', 'add');

    for j = 1:length(ind)
        if isfield(D,'Nfrequencies')
            h = imagesc(squeeze(D.data(handles.Cselection2(i), :,:, ind(j))));
            set(h, 'ButtonDownFcn', {@windowplot, i},...
                'Clipping', 'off');
        else
            h = plot([-D.events.start:D.events.stop]*1000/D.Radc,...
                D.data(handles.Cselection2(i), :, ind(j)),...
                'Color', handles.colour{j});
            set(h, 'ButtonDownFcn', {@windowplot, i},...
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
    handles.Heegfigures{ind} = [];
else
    handles.Heegfigures{ind} = figure;
    set(gcf,...
        'Name', (sprintf('Channel %s', D.channels.name{handles.Cselection2(ind)})),...
        'NumberTitle', 'off',...
        'DeleteFcn', {@delete_Heegwindows, ind});

    handles.Heegaxes2{ind} = gca;
    set(handles.Heegaxes2{ind}, 'NextPlot', 'add');

    if isfield(D,'Nfrequencies')
        xlabel('ms', 'FontSize', 16);
        ylabel('Hz', 'FontSize', 16, 'Interpreter', 'Tex')

        for i = 1:length(handles.Tselection)

            imagesc([-D.events.start:D.events.stop]*1000/D.Radc,D.tf.frequencies,squeeze(D.data(handles.Cselection2(ind), :,:, handles.Tselection(i))));

        end

        scale = get(handles.scaleslider, 'Value');
        set(gca, 'ZLim', [-scale scale],...
            'XLim',  [-D.events.start D.events.stop]*1000/D.Radc, 'YLim', [min(D.tf.frequencies) max(D.tf.frequencies)], 'Box', 'on');
        colormap('jet')
        caxis([-scale scale])

    else
        xlabel('ms', 'FontSize', 16);
        ylabel(D.units, 'FontSize', 16, 'Interpreter', 'Tex')

        for i = 1:length(handles.Tselection)
            plot([-D.events.start:D.events.stop]*1000/D.Radc,...
                D.data(handles.Cselection2(ind), :, handles.Tselection(i)),...
                'Color', handles.colour{i}, 'LineWidth', 2);
        end

        scale = get(handles.scaleslider, 'Value');

        set(gca, 'YLim', [-scale scale],...
            'XLim', [-D.events.start D.events.stop]*1000/D.Radc, 'Box', 'on');
        grid on
    end
    title(sprintf('%s (%d)', D.channels.name{handles.Cselection2(ind)}, handles.Cselection2(ind)), 'FontSize', 16);

    if D.Nevents > 1
        legend(handles.trialnames{handles.Tselection}, 0);
    end

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

% remove gfx struct
D = rmfield(D, 'gfx');

if spm_matlab_version_chk('7') >= 0
    save(fullfile(D.path, D.fname), '-V6', 'D');
else
    save(fullfile(D.path, D.fname), 'D');
end

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

if D.events.reject(ind) == 0
    % reject this trial
    tmp = [sprintf('%-12s', sprintf('trial %d', ind))  sprintf('%-4s', sprintf('%d', D.events.code(ind)))];
    tmp = [tmp sprintf('%-8s', 'reject')];
    D.events.reject(ind) = 1;
else
    % un-reject this trial
    tmp = [sprintf('%-12s', sprintf('trial %d', ind))  sprintf('%-4s', sprintf('%d', D.events.code(ind)))];
    D.events.reject(ind) = 0;
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

D.gfx.channels = ind;
S.D = D;
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
        if t(n) >= -D.events.start*1000/D.Radc...
                & t(n) <= D.events.stop*1000/D.Radc...
                & (strcmpi(s, '2d') | strcmpi(s, '3d'))
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
    % change time for james' function
    t=round(t/1000*D.Radc)+D.events.start+1;
    if length(t) == 1
        d = squeeze(D.data(D.channels.eeg, t, handles.Tselection(1)));
    else
        d = squeeze(mean(D.data(D.channels.eeg, t, handles.Tselection(1)), 2));
    end
    if strmatch(D.channels.ctf,'bdf_setup.mat')
        spm_eeg_scalp3d(d);
    else
        errordlg({'This can only be used for 128 channel BDF data at the moment'});
        return;
    end


end
spm('Pointer', 'Arrow');drawnow;

