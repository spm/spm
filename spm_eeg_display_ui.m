function Heeg = spm_eeg_display_ui(varargin)
% user interface for displaying EEG/MEG channel data.
% Heeg = spm_eeg_display_ui(varargin)
% 
% action    - first argument is a string that determines action taken
% optional 2nd argument:
% S		    - struct
%    fields of S:
%     D 	- EEG struct
%     Hfig	- Figure (or axes) to work in (Defaults to gcf)
%
% output:
%     Heeg		- Handle of resulting figure
%_______________________________________________________________________
% @(#)spm_eeg_display_ui.m	1.1 Stefan Kiebel 04/06/28

if nargin == 0
	error('Insufficient arguments');
elseif ~isstr(varargin{1})
	error('First argument must be action string');
end

action = varargin{1};

if nargin == 2
    S = varargin{2};
end

try
    D = S.D;
catch
    D = spm_get(1, '.mat', 'Select EEG mat file');
    try
        D = spm_eeg_ldata(D);
    catch    
        error(sprintf('Trouble reading file %s', D));
    end
end

if D.Nevents == 1 & ~isfield(D.events, 'start')
    errordlg({'Continuous data cannot be displayed.','Epoch first.'});
    return;
end

switch lower(action)
    
    %=======================================================================
    case {'setup', 'display'}
        %=======================================================================
        
        try
            F = S.Hfig;
        catch
            F = findobj('Tag', 'Graphics');
        end
               
end

% set(F, 'Tag', 'Graphics');

%         try
%             F = S.F;
%             
%             if isstr(F) 
%                 F = spm_figure('FindWin', F); 
%             end
%             
%             if ~ishandle(F)
%                 error('Invalid handle')
%             end
%             
%             switch get(F, 'Type')
%                 case 'figure'
%                     Heeg = [];
%                 case 'axes'
%                     Heegax = F;
%                     F = get(Heegax, 'Parent');
%                 otherwise
%                     error('F is not a figure or axis handle')
%             end
%         catch
%             [Finter, F, CmdLine] = spm('FnUIsetup','Display',0);
%             Heeg = [];
%         end

if strcmpi('setup', action)
    D.gfx = [];
    D.gfx.first = 1; % flag for first display
    D.gfx.Ncheckboxes = min(30, D.Nevents);
    D.gfx.linestyle = {'-', '--', '-.', ':'};
    D.gfx.colour = {[1 0 0], [0 1 0], [0 0 1], [1 1 0]};
    D.gfx.showCurrent = 1;
    
    D.gfx.trial = 1;
    
    if D.gfx.showCurrent == 1
        D.gfx.Tdisplay = D.gfx.trial; 
    else
        D.gfx.Tdisplay = []; 
    end
    
    D.gfx.Cdisplay = []; % display of channels
    
    D.gfx.Rcolour = {[0 1 0], [1 0 0]}; % colour coding for accept, reject
    
    set(F, 'UserData', []);
    figure(F);clf
    
    % text in first line
    % 		uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
    % 			'Position', [0.02 0.95 0.5 0.03],...
    % 			'String', ['Time stamp ', datestr([D.descrip.date, ' ', D.descrip.time], 0)],...
    % 			'BackgroundColor', get(F, 'Color'), 'HorizontalAlignment', 'left');
    
    % Slider for trial number
    FS1 = spm('FontSize', 16);
    if D.Nevents > 1
        h = uicontrol(F, 'Tag', 'trialselect', 'Style', 'slider',...
            'Min', 1, 'Max', D.Nevents, 'Value', 1, 'Units',...
            'normalized', 'Position', [0.05 0.36 0.25 0.03],...
            'SliderStep', [1/(D.Nevents-1) min(D.Nevents-1, 10/(D.Nevents-1))],...
            'TooltipString', 'Choose trial number',...
            'Callback', 'spm_eeg_update_gfx(''trial'', round(get(gcbo, ''Value'')))',...
            'Parent', F, 'Interruptible', 'off');
    uicontrol(F, 'Style','Frame','BackgroundColor',spm('Colour'), 'Units',...
        'normalized', 'Position',[0.05 0.32 0.25 0.04]);
    
    uicontrol(F, 'Style', 'text', 'Units', 'normalized', 'String', '1',...
        'Position',[0.06 0.325 0.06 0.03],...
        'HorizontalAlignment', 'left', 'FontSize', FS1,...
        'BackgroundColor',spm('Colour'));
    
    uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
        'String', mat2str(D.Nevents),...
        'Position',[0.23 0.325 0.06 0.03],...
        'HorizontalAlignment', 'right', 'FontSize', FS1,...
        'BackgroundColor',spm('Colour'));
    
    uicontrol(F, 'Style', 'text', 'Tag', 'trialtext',...
        'Units', 'normalized',...
        'String', mat2str(round(get(findobj('Tag', 'trialselect'), 'Value'))),...
        'Position',[0.14 0.325 0.06 0.03],...
        'HorizontalAlignment', 'center', 'FontSize', FS1,...
        'BackgroundColor',spm('Colour'));
            end

    % Compute first scaling value
    D.gfx.scale = 4*ceil(max(max(abs(D.data(:, :, 1)))));
    
    uicontrol('Tag', 'scaleselect', 'Style', 'slider',...
        'Min', 1, 'Max', D.gfx.scale, 'Value', D.gfx.scale/2, 'Units',...
        'normalized', 'Position', [0.4 0.36 0.25 0.03],...
        'SliderStep', [1/(D.gfx.scale-1) 10/(D.gfx.scale-1)],...
        'TooltipString', 'Choose scaling',...
        'Callback', 'spm_eeg_update_gfx(''scaling'', get(gcbo, ''Value''))',...
        'Parent', F);

    uicontrol(F, 'Style','Frame','BackgroundColor',spm('Colour'), 'Units',...
        'normalized', 'Position',[0.40 0.32 0.25 0.04]);
    
    uicontrol(F, 'Style', 'text', 'Units', 'normalized', 'String', '1',...
        'Position',[0.41 0.325 0.11 0.03],...
        'HorizontalAlignment', 'left', 'FontSize', FS1,...
        'BackgroundColor',spm('Colour'));
    
    uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
        'String', mat2str(D.gfx.scale),...
        'Position',[0.53 0.325 0.11 0.03],...
        'HorizontalAlignment', 'right', 'FontSize', FS1,...
        'BackgroundColor',spm('Colour'));
    
    uicontrol(F, 'Style', 'text', 'Tag', 'scaletext',...
        'Units', 'normalized',...
        'String', mat2str(D.gfx.scale/2),...
        'Position',[0.49 0.325 0.06 0.03],...
        'HorizontalAlignment', 'center', 'FontSize', FS1,...
        'BackgroundColor',spm('Colour'));
    
    % Setup togglebutton for accept/reject
    uicontrol('Tag', 'rejectbutton', 'Style', 'togglebutton',...
        'Units', 'normalized', 'Position', [0.7 0.4 0.1 0.05],...
        'String', 'Reject', 'FontSize', 16, 'Min', 0, 'Max', 1,...
        'SelectionHighlight', 'off', 'BackgroundColor', D.gfx.Rcolour{D.events.reject(D.gfx.trial)+1},...
        'CallBack', 'spm_eeg_update_gfx(''reject'', get(gcbo, ''Value''))',...
        'Parent', F);
    
    % Setup pushbutton for save
    uicontrol('Tag', 'savebutton', 'Style', 'pushbutton',...
        'Units', 'normalized', 'Position', [0.7 0.3 0.1 0.05],...
        'String', 'Save', 'FontSize', 16,...
        'BackgroundColor', spm('Colour'),...
        'CallBack', 'spm_eeg_update_gfx(''save'')',...
        'Parent', F);
    
    % setup trial display slider
    if D.Nevents > 1
        uicontrol('Tag', 'rejectslider', 'Style', 'slider',...
            'Units', 'normalized', 'Position', [0.94 0.3 0.03 0.02*D.gfx.Ncheckboxes],...
            'Min', 1, 'Max', D.Nevents,...
            'SliderStep', [1/(D.Nevents-1) min(D.Nevents-1, 10/(D.Nevents-1))],...
            'Value', get(findobj('Tag', 'trialselect'), 'Value'),...
            'CallBack', 'spm_eeg_update_gfx(''trial'', round(get(gcbo, ''Value'')))',...
            'Parent', F, 'Interruptible', 'off');
        
        % setup trial pushbutton list
        for i = 1:D.gfx.Ncheckboxes
            h = uicontrol('Tag', sprintf('displaylist%d', i), 'Style', 'checkbox',...
                'Units', 'normalized', 'Position', [0.85 0.28+i*0.02 0.09 0.02],...
                'Min', 0, 'Max', 1, 'String', mat2str(i),...
                'HorizontalAlignment', 'left',...
                'Value', ismember(i, D.gfx.Tdisplay),...
                'BackgroundColor', D.gfx.Rcolour{D.events.reject(i)+1},...
                'CallBack', 'spm_eeg_update_gfx(''display'', get(gcbo, ''Value''))',...
                'Parent', F);
            
            if i == get(findobj('Tag', 'trialselect'), 'Value')
                set(h, 'FontWeight', 'bold');
            end
        end
    end
    cmenu = uicontextmenu('Parent', F, 'Tag', 'bigplotmenu');
    uimenu(cmenu, 'label', 'Clear graph', 'CallBack', 'spm_eeg_update_gfx(''cleargraph'')');
    
    % setup larger plot of selected time-series at the bottom
    h = axes('Tag', 'bigplot', 'Position', [0.15 0.08 0.65 0.2], 'Visible', 'off',...
        'Parent', F, 'XLim', [-D.events.start D.events.stop],...
        'YLim', [-D.gfx.scale D.gfx.scale],...
        'UIContextMenu', cmenu);
    
    
end

% Display flathead image
%-----------------------------------------------------------------------
clear S
S.D = D;
S.Hfig = F;
% set(F, 'UserData', D);

D = spm_eeg_display(S);



Heeg = {F};

