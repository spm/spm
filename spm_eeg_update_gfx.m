function spm_eeg_update_gfx(varargin)
% (internal) function to display epoched EEG/MEG channel data
% FORMAT spm_eeg_update_gfx(S)
%
% action    - first argument is a string that determines action taken
%_______________________________________________________________________
%
% spm_eeg_update_gfx is an internally used function that plots EEG/MEG
% traces.
%_______________________________________________________________________
% %W% Stefan Kiebel %E%

% get Userdata
F = findobj('Tag', 'Graphics');
D = get(F, 'UserData');

switch(lower(varargin{1}))
    
    %=======================================================================
    case 'reject'
        %=======================================================================
        s = varargin{2};
        t = D.gfx.trial;
        
        set(gcbo, 'BackgroundColor', D.gfx.Rcolour{s+1});
        
        set(findobj('Style', 'checkbox', 'String', num2str(t)), 'BackgroundColor', D.gfx.Rcolour{s+1});		
        
        D.events.reject(t) = s;
        
        set(F, 'UserData', D);
        
        %=======================================================================
    case 'display'
        %=======================================================================
        % called by checkbuttons only...
        s = varargin{2};
        
        % which trial?
        t = str2num(get(gcbo, 'String'));
        
        if D.gfx.showCurrent == 1 & t == D.gfx.trial
            set(gcbo, 'Value', 1);
        else
            if length(D.gfx.Tdisplay) < length(D.gfx.linestyle)
                if s == 0
                    D.gfx.Tdisplay(D.gfx.Tdisplay == t) = [];
                    D.gfx.NextPlot = 'replace';
                else
                    D.gfx.Tdisplay = unique([D.gfx.Tdisplay t]);
                    D.gfx.NextPlot = 'replace';
                end
            else
                set(gcbo, 'Value', 1);
            end
            
            clear S; S.D = D; S.Hfig = F;
            spm_eeg_display(S);
            
            if ~isempty(D.gfx.Cdisplay) & ~isempty(D.gfx.Tdisplay)
                makebigplot(D);
            else 
                spm_eeg_update_gfx('cleargraph');
            end
        end
        
        %=======================================================================
    case 'trial'
        %=======================================================================
        
        t = varargin{2};
        
        if D.gfx.showCurrent == 1
            % remove old trial from display list
            D.gfx.Tdisplay(D.gfx.Tdisplay == D.gfx.trial) = [];
            % add new trial to display list
            D.gfx.Tdisplay = unique([D.gfx.Tdisplay t]);
        end
        
        D.gfx.NextPlot = 'replace';
        
        D.gfx.trial = t;
        
        % Update trial list, keep (if possible) current trial in centre
        Ncb = D.gfx.Ncheckboxes;
        c = floor(Ncb/2);
        
        for i = 1:Ncb
            if t <= c
                set(findobj('Tag', sprintf('displaylist%d', i)), 'String', mat2str(i),...
                    'Value', ismember(i, D.gfx.Tdisplay),...
                    'BackgroundColor', D.gfx.Rcolour{D.events.reject(i)+1},...
                    'FontWeight', 'normal');			
                if i == t
                    set(findobj('Tag', sprintf('displaylist%d', i)), 'FontWeight', 'bold');
                end
            elseif t > D.Nevents - Ncb + c
                set(findobj('Tag', sprintf('displaylist%d', i)), 'String', mat2str(D.Nevents-Ncb+i),...
                    'Value', ismember(D.Nevents-Ncb+i, D.gfx.Tdisplay),...
                    'BackgroundColor', D.gfx.Rcolour{D.events.reject(D.Nevents-Ncb+i)+1},...
                    'FontWeight', 'normal');
                if i == t - D.Nevents + Ncb
                    set(findobj('Tag', sprintf('displaylist%d', i)), 'FontWeight', 'bold');
                end
            else
                set(findobj('Tag', sprintf('displaylist%d', i)), 'String', mat2str(t-c+i),...
                    'Value', ismember(t-c+i, D.gfx.Tdisplay),...
                    'BackgroundColor', D.gfx.Rcolour{D.events.reject(t-c+i)+1},...
                    'FontWeight', 'normal');
                if i == c
                    set(findobj('Tag', sprintf('displaylist%d', i)), 'FontWeight', 'bold');
                end
            end
        end
        
        drawnow
        
        % Update various displays of the trial number
        set(findobj(F, 'Tag', 'trialtext'), 'String', mat2str(t));
        
        set(findobj(F, 'Tag', 'trialselect'), 'Value', t);
        
        set(findobj(F, 'Tag', 'rejectslider'), 'Value', t);
        
        
        % Update colour of reject button
        if D.events.reject(t) == 0
            set(findobj(F, 'Tag', 'rejectbutton'),...
                'Value', 0, 'BackgroundColor', [0 1 0]);
        else
            set(findobj(F, 'Tag', 'rejectbutton'),...
                'Value', 1, 'BackgroundColor', [1 0 0]);
        end
        
        % Update main display
        clear S; S.D = D; S.Hfig = F;
        spm_eeg_display(S);
      
        if ~isempty(D.gfx.Cdisplay)
            makebigplot(D);
        end
        %=======================================================================
    case 'scaling'
        %=======================================================================
        s = varargin{2};
        D.gfx.scale = s;
        
        clear S; S.D = D; S.F = F;
        spm_eeg_display_ui('display', S);
        
        set(findobj(F, 'Tag', 'scaletext'), 'String', mat2str(round(D.gfx.scale)));
        
        set(findobj(F, 'Tag', 'scaletext2'), 'String', sprintf(' %d \\muV', round(D.gfx.scale)));
        
        if ~isempty(D.gfx.Cdisplay)
            makebigplot(D);
        end
        %=======================================================================
    case 'addchannel'
        %=======================================================================
        % adds channeldata to big graph
        h = varargin{2}; % handle of uimenu
        c = get(h, 'UserData'); % channel to add
        
        if ~ismember(c, D.gfx.Cdisplay) & length(D.gfx.Cdisplay) < length(D.gfx.colour)
            D.gfx.Cdisplay = unique([D.gfx.Cdisplay c]);
            makebigplot(D);	
            
            
            % Change contextmenu of bigplot
            uimenu(findobj('Tag', 'bigplotmenu'), 'label',...
                sprintf('Remove %s', D.channels.name{c}), 'CallBack', 'spm_eeg_update_gfx(''rmchannel'', gcbo)',...
                'UserData', c);
            
            set(findobj('Tag', 'bigplot'), 'Visible', 'on');
            
            % Update menu of small plot
            set(D.gfx.Heegmenus_rm(c), 'Enable', 'on');
            set(D.gfx.Heegmenus_add(c), 'Enable', 'off');
            
            set(F, 'UserData', D);
        end
        %=======================================================================
    case 'rmchannel'
        %=======================================================================
        % remove channeldata from big graph
        h = varargin{2}; % handle of uimenu
        c = get(h, 'UserData'); % channel to remove
        
        D.gfx.Cdisplay = setxor(D.gfx.Cdisplay, c);
        
        if isempty(D.gfx.Cdisplay)
            spm_eeg_update_gfx('cleargraph');
        else
            makebigplot(D);
        end
        
        % remove uimenu from uicontextmenu
        h = findobj('Label', sprintf('Remove %s', D.channels.name{c}));
        delete(h);
        
        % Update menu of small plot
        set(D.gfx.Heegmenus_rm(c), 'Enable', 'off');
        set(D.gfx.Heegmenus_add(c), 'Enable', 'on');
        
        set(F, 'UserData', D);
        %=======================================================================
    case 'cleargraph'
        %=======================================================================
        D.gfx.Cdisplay = [];
        set(D.gfx.Heegmenus_rm, 'Enable', 'off');
        set(D.gfx.Heegmenus_add, 'Enable', 'on');
        
        
        h = findobj('Tag', 'bigplot');
        
        % Kill all children (except for the last) of uicontrolmenu
        c = get(findobj('Tag', 'bigplotmenu'), 'Children');
        
        delete(c(1:end-1));
        
        axes(h);
        cla;
        legend('off');
        set(h, 'Visible', 'off');
        set(F, 'UserData', D);

        %=======================================================================
    case 'save'
        %=======================================================================
        spm('Pointer', 'Watch');
        gfx = D.gfx;
        D = rmfield(D, 'gfx');
        
        save(fullfile(D.path, D.fname), 'D');
        D.gfx = gfx;
        spm('Pointer', 'Arrow');
        
    otherwise
        error('Unknown option: %s', varargin{1});
        
end % switch


function makebigplot(D)

Ci = length(D.gfx.Cdisplay);
h = findobj('Tag', 'bigplot');
axes(h);
cla
set(h, 'NextPlot', 'add');
xlabel('ms', 'FontSize', 16); ylabel('\muV', 'FontSize', 16, 'Interpreter', 'Tex')

for i = 1:length(D.gfx.Tdisplay)
    for j = 1:length(D.gfx.Cdisplay)
        % plot data as function of peristimulus time
        plot([-D.events.start:D.events.stop]*1000/D.Radc, D.data(D.gfx.Cdisplay(j), :, D.gfx.Tdisplay(i)), 'Color', D.gfx.colour{j}, 'LineStyle', D.gfx.linestyle{i});
    end
end
set(gca, 'YLim', [-D.gfx.scale/2 D.gfx.scale/2],...
    'XLim', [-D.events.start D.events.stop]*1000/D.Radc, 'Box', 'on');
grid on

legend(D.channels.name{D.gfx.Cdisplay}, 0);

