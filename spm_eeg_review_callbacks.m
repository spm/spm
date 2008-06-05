function [] = spm_eeg_review_callbacks(arg1,arg2,arg3)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jean Daunizeau
% $Id: spm_eeg_review_callbacks.m 1792 2008-06-05 15:43:54Z jean $

D = get(gcf,'userdata');
handles = D.PSD.handles;

switch arg1

    
    %% File I/O 
    case 'file'
        switch arg2
            case 'save'
                D = rmfield(D,'PSD');
                save(D.fname,'D')
        end
        
        
    
    
    %% Visualization callbacks

    case 'visu'

        zoom off

        switch arg2

            case 'main'

                try
                    handles = rmfield(D.PSD.handles,'PLOT');
                    D.PSD.handles = handles;
                end
                switch arg3
                    case 1
                        [D] = spm_eeg_review_switchDisplay(D,'standardData');
                        if strcmp(D.PSD.VIZU.type,'standardData')
                            updateDisp(D)
                        end
                    case 2
                        [D] = spm_eeg_review_switchDisplay(D,'scalpData');
                        if strcmp(D.PSD.VIZU.type,'scalpData')
                            D.PSD.trials.current = D.PSD.trials.current(1);
                            updateDisp(D)
                        end
                    case 3
                        [D] = spm_eeg_review_switchDisplay(D,'visuRecon');
                        if strcmp(D.PSD.VIZU.type,'visuRecon')
                            updateDisp(D)
                        end
                end
                

                
            case 'update'
                
                updateDisp(D)
                
                % Scalp interpolation
            case 'scalp_interp'

                if ~isempty([D.channels(:).X_plot2D])
                    x = round(mean(get(handles.axes(1),'xlim')));
                    ylim = get(handles.axes,'ylim');
                    if strcmp(D.PSD.VIZU.type,'standardData')
                        hl = line('parent',handles.axes,'xdata',[x;x],'ydata',[ylim(1);ylim(2)]);
                        in.hl = hl;
                    end
                    switch D.PSD.type
                        case 'continuous'
                            trN = 1;
                            in.gridTime = [1:D.Nsamples]./D.Fsample + D.timeOnset;
                            in.unit = 's';
                        case 'epoched'
                            trN = D.PSD.trials.current(1);
                            in.trN = trN;
                            in.gridTime = [1:D.Nsamples].*1e3./D.Fsample + D.timeOnset.*1e3;
                            in.unit = 'ms';
                    end
                    in.x = x;

                    in.handles = handles;
                    % for EEG and MEG sensors
                    Ieeg  = find(strcmp('EEG',{D.channels.type}));
                    Ieeg = intersect(Ieeg,find(~[D.channels.bad]));
                    if ~isempty(Ieeg)
                        posEEG(:,1) = [D.channels(Ieeg).X_plot2D]';
                        posEEG(:,2) = [D.channels(Ieeg).Y_plot2D]';
                        EEGlabels = {D.channels(Ieeg).label};
                        yEEG = D.data.y(Ieeg,:,trN);
                        in.min = min(yEEG(:));
                        in.max = max(yEEG(:));
                        in.ind = Ieeg;
                        in.type = 'EEG';
                        yEEG = yEEG(:,x);
                        spm_eeg_plotScalpData(yEEG,posEEG,EEGlabels,in);
                    end
                    Imeg  = find(strcmp('MEG',{D.channels.type}));
                    Imeg = intersect(Imeg,find(~[D.channels.bad]));
                    if ~isempty(Imeg)
                        posMEG(:,1) = [D.channels(Imeg).X_plot2D]';
                        posMEG(:,2) = [D.channels(Imeg).Y_plot2D]';
                        MEGlabels = {D.channels(Imeg).label};
                        yMEG = D.data.y(Imeg,:,trN);
                        in.min = min(yMEG(:));
                        in.max = max(yMEG(:));
                        in.ind = Imeg;
                        in.type = 'MEG';
                        yMEG = yMEG(:,x);
                        spm_eeg_plotScalpData(yMEG,posMEG,MEGlabels,in);
                    end
                else
                    msgbox('Get 2d positions for EEG/MEG channels!')
                end

            case 'inv'
                

                if ~isequal(D.PSD.invN,arg3)
                    delete(findobj('tag','dipSpheres'))
                    str = getInfo4Inv(D,arg3);
                    isInv = get(D.PSD.handles.BMCcurrent,'userdata');
                    set(D.PSD.handles.infoText,'string',str);
                    set(D.PSD.handles.BMCcurrent,'XData',find(isInv==arg3));
                    D.PSD.invN = arg3;
                    trN = D.PSD.trials.current;
                    model = D.other.inv{D.PSD.invN}.inverse;
                    J = model.J{trN}*model.T';
                    set(D.PSD.handles.axes,'CLim',[min(min(J)) max(max(J))]);
                    set(D.PSD.handles.mesh,...
                        'Vertices',D.other.inv{D.PSD.invN}.mesh.tess_mni.vert,...
                        'Faces',D.other.inv{D.PSD.invN}.mesh.tess_mni.face);
                    if isfield(D.other.inv{D.PSD.invN}.inverse,'dipfit')
                        xyz = D.other.inv{D.PSD.invN}.inverse.dipfit.Lpos;
                        Np  = size(xyz,2);
                        radius = D.other.inv{D.PSD.invN}.inverse.dipfit.radius;
                        [x,y,z] = sphere(20);
                        axes(D.PSD.handles.axes)
                        for i=1:Np
                            D.PSD.handles.dipSpheres(i) = patch(...
                                surf2patch(x.*radius+xyz(1,i),...
                                y.*radius+xyz(2,i),z.*radius+xyz(3,i)));
                            set(D.PSD.handles.dipSpheres(i),'facecolor',[1 1 1],...
                                'edgecolor','none','facealpha',0.5,...
                                'tag','dipSpheres');
                        end
                    end
                    
                    updateDisp(D);
                end


                % Contrast/intensity rescaling
            case 'iten_sc'

                D.PSD.VIZU.visu_scale = arg3*D.PSD.VIZU.visu_scale;
                if strcmp(D.PSD.VIZU.type,'standardData')
                    D.PSD.VIZU.xlim = get(D.PSD.handles.axes(1),'xlim');
                end
                updateDisp(D);




                % Resize plotted data window
            case 'time_w'

                % Get current plotted data window range and limits
                xlim = get(handles.axes(1),'xlim');
                %                 ylim = get(handles.axes(1),'ylim');
                length_window = max([arg3*round(xlim(2)-xlim(1)),1]);
                xm = mean(xlim);

                % Change limits of plotted data window
                xlim = round([xm-length_window./2 , xm+length_window./2]);
                xlim(1) = max([xlim(1) 1]);
                if length_window >= D.Nsamples-1
                    xlim = [1 D.Nsamples];
                else
                    if xlim(2) >=D.Nsamples
                        dx = xlim(2) - D.Nsamples;
                        xlim(2) = D.Nsamples;
                        xlim(1) = xlim(1) - dx -1;
                        xlim(1) = max([xlim(1) 1]);
                    end
                end
                D.PSD.VIZU.xlim = xlim;

                % This part avoids limiting displaying conditions
                if length_window >= D.Nsamples-1
%                     set(handles.VIZU.time_w1,'enable','off')
                    set(handles.BUTTONS.vb3,'enable','off')
                elseif length_window < 20
%                     set(handles.VIZU.time_w2,'enable','off')
                    set(handles.BUTTONS.vb4,'enable','off')
                end

                % This part fixes the boundaries exceptions
                length_window = xlim(2)-xlim(1);
                ratio = length_window/200;
                maxi = D.Nsamples-length_window+1;
                val = get(handles.BUTTONS.slider_step,'value');
                if val > maxi;
                    set(handles.BUTTONS.slider_step,'value',maxi);
                end
                if arg3 > 1
                    set(handles.BUTTONS.slider_step,'visible','on');
                end
                if isequal(xlim,[1 D.Nsamples]) == 0
                    set(handles.BUTTONS.slider_step,...
                        'sliderstep',[ratio*10/(D.Nsamples-1) ratio*20/(D.Nsamples-1)],...
                        'visible','on');
                    set(D.PSD.handles.BUTTONS.goPlusOne,'visible','on');
                    set(D.PSD.handles.BUTTONS.goMinusOne,'visible','on');
                else
                    set(handles.BUTTONS.slider_step,'visible','off');
                    set(D.PSD.handles.BUTTONS.goPlusOne,'visible','off');
                    set(D.PSD.handles.BUTTONS.goMinusOne,'visible','off');
                end

                if arg3 > 1
%                     set(handles.VIZU.time_w2,'enable','on');
                    set(handles.BUTTONS.vb4,'enable','on');
                else
%                     set(handles.VIZU.time_w1,'enable','on');
                    set(handles.BUTTONS.vb3,'enable','on');
                end

                updateDisp(D,1)



                % Zoom (box in)
            case 'zoom'
                
                switch D.PSD.VIZU.type

                    case 'standardData'

                        if arg3
                            zoom
                            %                     set(handles.VIZU.zoom2,'enable','on')
                        else % reset zoom and rebuild normal plotted data window
                            %                     set(handles.VIZU.zoom2,'enable','off')
                            updateDisp(D)
                        end

                    case 'scalpData'
                        
                        try
                            axes(D.PSD.handles.scale)
                        end
                        [x, y, button] = ginput(1);
                        indAxes = get(gco,'userdata');
                        if ~~indAxes
                            hf = figure;
                            chanLabel = D.channels(D.PSD.VIZU.visuSensors(indAxes)).label;
                            if D.channels(D.PSD.VIZU.visuSensors(indAxes)).bad
                                chanLabel = [chanLabel,' (BAD)'];
                            end
                            set(hf,'name',['channel ',chanLabel])
                            ha2 = axes('parent',hf,...
                                'xlim',get(D.PSD.handles.axes(indAxes),'xlim'),...
                                'ylim',get(D.PSD.handles.axes(indAxes),'ylim'),...
                                'nextplot','add',...
                                'XGrid','on','YGrid','on');
                            trN = D.PSD.trials.current(:);
                            Ntrials = length(trN);
                            
                            if strcmp(D.transform.ID,'time')

                                leg = cell(Ntrials,1);
                                col = colormap('lines');
                                col = repmat(col(1:7,:),floor(Ntrials./7)+1,1);
                                hp = get(handles.axes(indAxes),'children');
                                for i=1:Ntrials
                                    datai = get(hp(Ntrials-i+1),'ydata')./D.PSD.VIZU.visu_scale;
                                    hp2(i) = plot(ha2,datai,...
                                        'color',col(i,:));
                                    leg{i} = D.PSD.trials.TrLabels{trN(i)};
                                end
                                legend(leg)
                                xg = 0:D.Fsample/10:D.Nsamples;
                                set(gca,'xtick',xg,'xticklabel',xg./D.Fsample*1e3+D.timeOnset*1e3);
                                xlabel(ha2,'time (in ms after time onset)')
                                title(ha2,['channel ',chanLabel,...
                                    ' (',D.channels(D.PSD.VIZU.visuSensors(indAxes)).type,')'])

                            else

                                miY = 0;
                                maY = 0;

                                datai = squeeze(D.data.y(indAxes,:,:,trN(1)));
                                hp2 = imagesc(datai);
                                set(hp2,'parent',ha2);
                                colormap('jet')
                                colorbar
                                xg = 0:D.Fsample/10:D.Nsamples;
                                set(gca,'xtick',xg,'xticklabel',xg./D.Fsample*1e3+D.timeOnset*1e3);
                                xlabel(ha2,'time (in ms after time onset)')
                                ytick = get(ha2,'ytick');
                                set(ha2,'yticklabel',D.transform.frequencies(ytick))
                                ylabel(ha2,'frequency (in Hz)')
                                title(ha2,['channel ',chanLabel,...
                                    ' (',D.channels(D.PSD.VIZU.visuSensors(indAxes)).type,')'])
                                
                            end
                            
                            axes(ha2)
                        end
                end


                % Select all sensors
            case 'sensor_select_all'

                D.PSD.VIZU.visuSensors = 1:length(D.channels);
                updateDisp(D)


                %% select sensors using dedicated GUI
            case 'sensor_select'

                in = PSD_gui_selectSensors(D,1);
                if ~isequal(in,1:length(D.channels))
                    % Get select sensors
                    D.PSD.VIZU.visuSensors = in;
                    updateDisp(D)
                end



%                 %% Data navigation using the editable box
%             case 'focus_t'
% 
% 
%                 % Get current plotted data window range and limits
%                 xlim0 = get(handles.axes,'xlim');
%                 
%                 val = str2num(...
%                     get(D.PSD.handles.BUTTONS.slider_step,'string'));
% 
%                 offset = round(val);
% 
%                 try
%                     % The IF statement ensures acceptable range
%                     if ~isequal(xlim0,[1 D.Nsamples])
% 
%                         % Build limits of the plotted data window
%                         length_window = round(xlim0(2)-xlim0(1));
%                         if offset < round(0.5*length_window)
%                             offset = round(0.5*length_window);
%                             set(handles.BUTTONS.slider_step,'value',1);
%                         elseif offset > D.Nsamples-round(0.5*length_window)
%                             offset = D.Nsamples-round(0.5*length_window)-1;
%                             set(handles.BUTTONS.slider_step,'value',get(handles.BUTTONS.slider_step,'max'));
%                         else
%                             set(handles.BUTTONS.slider_step,'value',offset);
%                         end
%                         xlim = [offset-round(0.5*length_window) offset+round(0.5*length_window)];
%                         xlim(1) = max([xlim(1) 1]);
%                         xlim(2) = min([xlim(2) D.Nsamples]);
%                         set(gco,'string',num2str(mean(xlim)));
% 
%                         D.PSD.VIZU.xlim = xlim;
%                         D.PSD.VIZU.x0 = offset;
%                         
%                         updateDisp(D,1)
% 
%                     end
%                 catch
%                     set(D.PSD.handles.BUTTONS.slider_step,'string',num2str(mean(xlim0)));
%                     spm_eeg_review_callbacks('visu','focus_t',0);
%                 end


                %% Data navigation using the slider
            case 'slider_t'
                
                offset = get(gco,'value');
                if ~strcmp(D.PSD.VIZU.type,'visuRecon')
                    offset = round(offset);
                    % Get current plotted data window range and limits
                    xlim0 = get(handles.axes(1),'xlim');
                    % The IF statement ensures acceptable range
                    if isequal(xlim0,[1 D.Nsamples]) == 0
                        % deal w/ boundaries of the dataset window
                        length_window = xlim0(2)-xlim0(1);
                        xlim = round([offset-length_window/2 offset+length_window/2]);
                        if xlim(1) < 1 && xlim(2) <= D.Nsamples
                            xlim(1) = 1;
                            xlim(2) = round(1 + length_window);
                            set(handles.BUTTONS.slider_step,'value',mean(xlim));
                        elseif xlim(1) >= 1 && xlim(2) > D.Nsamples
                            xlim(2) = D.Nsamples;
                            xlim(1) = round(D.Nsamples - length_window);
                            set(handles.BUTTONS.slider_step,'value',mean(xlim));
                        elseif xlim(1) < 1 && xlim(2) > D.Nsamples
                            xlim(1) = 1;
                            xlim(2) = D.Nsamples;
                            set(handles.BUTTONS.slider_step,'value',mean(xlim));
                        end
%                         set(handles.BUTTONS.focus_temp,'string',round(mean(xlim)));
                        D.PSD.VIZU.xlim = xlim;
                        updateDisp(D,1)
                    end
                else
                    D.PSD.VIZU.x0 = offset;
                    
                    updateDisp(D,1)
                end

                %% Scroll page by page (button)
            case 'goOne'

                % Get current plotted data window range and limits
                xlim0 = get(handles.axes(1),'xlim');
                xm = mean(xlim0);
                length_window = abs(diff(xlim0));
                if arg3 == 0
                    offset = xm - length_window;
                else
                    offset = xm + length_window;
                end
                % The IF statement ensures acceptable range
                if isequal(xlim0,[1 D.Nsamples]) == 0

                    % deal w/ boundaries of the dataset window
                    xlim = round([offset-length_window/2 offset+length_window/2]);
                    if xlim(1) < 1 && xlim(2) <= D.Nsamples
                        xlim(1) = 1;
                        xlim(2) = round(1 + length_window);
                    elseif xlim(1) >= 1 && xlim(2) > D.Nsamples
                        xlim(2) = D.Nsamples;
                        xlim(1) = round(D.Nsamples - length_window);
                    elseif xlim(1) < 1 && xlim(2) > D.Nsamples
                        xlim(1) = 1;
                        xlim(2) = D.Nsamples;
                    end
                    set(handles.BUTTONS.slider_step,'value',mean(xlim));
%                     set(handles.BUTTONS.focus_temp,'string',round(mean(xlim)));
                    D.PSD.VIZU.xlim = xlim;

                    updateDisp(D)

                end





                %% X/Y Grids
            case 'ygrid'

                gr = get(gca,'ygrid');
                if isequal(gr,'on')
                    set(gca,'ygrid','off');
                elseif isequal(gr,'off')
                    set(gca,'ygrid','on');
                end


            case 'xgrid'

                if isfield(D,'Fsample') && ~isempty(D.Fsample) && ~isequal(D.Fsample,0)
                    lab = get(gcbo,'label');
                    if isequal(lab,'x-axis grid: #seconds')
                        xg = 0:D.Fsample:D.Nsamples;
                        set(gca,'xtick',xg)
                        set(gca,'xgrid','on')
                        set(gcbo,'label','x-axis grid: #time samples')
                    elseif isequal(lab,'x-axis grid: #time samples')
                        xg = 0:500:D.Nsamples;
                        set(gca,'xtick',xg)
                        set(gca,'xgrid','on')
                        set(gcbo,'label','x-axis grid: #seconds')
                    end
                else
                    msgbox('No sampling frequency specified!');
                end



                %% Reverse data sign
            case  'MainSwitch'

                D.PSD.VIZU.montage.M = -D.PSD.VIZU.montage.M;
                updateDisp(D)


                %% other ?
            otherwise;disp('unknown command !')


        end


    case 'menuEvent'

        Nevents = length(D.trials.events);

        x                       = [D.trials.events.time]';
        x(:,2)                  = [D.trials.events.duration]';
        x(:,2)                  = sum(x,2);

        % Find the index of the selected event
        currentEvent = get(gco,'userdata');
        eventType = D.trials.events(currentEvent).type;
        eventValue = D.trials.events(currentEvent).value;
        tit = ['Current event is selection #',num2str(currentEvent),...
            ' /',num2str(Nevents),' (type= ',eventType,', value=',num2str(eventValue),').'];

        D.PSD.tools.undo.select = D.trials.events;


        switch arg2

            % Execute actions accessible from the event contextmenu : click
            case 'click'

                % Highlight the selected event
                hh = findobj('selected','on');
                set(hh,'selected','off');
                set(gco,'selected','on')

                % Prompt basic information on the selected event
                disp(tit)

                % Execute actions accessible from the event contextmenu : edit event properties
            case 'EventProperties'

                set(gco,'selected','on')

                % Build GUI for manipulating the event properties
                stc{1} = 'Current event is a selection of type...';
                stc{2} = 'Current event has value...';
                stc{3} = 'Starts at (sec)...';
                stc{4} = 'Duration (sec)...';
                default{1} = eventType;
                default{2} = num2str(eventValue);
                default{3} = num2str(x(currentEvent,1));
                default{4} = num2str(abs(diff(x(currentEvent,:))));
                answer = inputdlg(stc,tit,1,default);

                if ~isempty(answer)

                    try
                        eventType = answer{1};
                        eventValue = str2num(answer{2});
                        D.trials.events(currentEvent).time = str2num(answer{3});
                        D.trials.events(currentEvent).duration = str2num(answer{4});
                        D.trials.events(currentEvent).type = eventType;
                        D.trials.events(currentEvent).value = eventValue;
                    end

                    D.PSD.tools.redo.select = D.trials.events;

                    handles = rmfield(D.PSD.handles,'PLOT');
                    D.PSD.handles = handles;
                    updateDisp(D)

                    %                     set(handles.EDIT.undo,'enable','on');
                    %                     set(handles.EDIT.redo,'enable','off');

                end


                % Execute actions accessible from the event contextmenu : go to next/previous event
            case 'goto'


                here = mean(x(currentEvent,:));

                values = [D.trials.events.value];
%                 sameValue = find(values==eventValue);
%                 xm = mean(x(sameValue,:),2);
                xm = mean(x(values==eventValue,:),2);
                if arg3 == 0
                    ind = find(xm < here);
                else
                    ind = find(xm > here);
                end

                if ~isempty(ind)

                    if arg3 == 0
                        offset = round(max(xm(ind))).*D.Fsample;
                    else
                        offset = round(min(xm(ind))).*D.Fsample;
                    end
%                     ud = get(handles.axes,'userdata');

%                     ylim = D.PSD.VIZU.ylim;
                    xlim0 = get(handles.axes,'xlim');

                    if ~isequal(xlim0,[1 D.Nsamples])

                        length_window = round(xlim0(2)-xlim0(1));
                        if offset < round(0.5*length_window)
                            offset = round(0.5*length_window);
                            set(handles.BUTTONS.slider_step,'value',1);
                        elseif offset > D.Nsamples-round(0.5*length_window)
                            offset = D.Nsamples-round(0.5*length_window)-1;
                            set(handles.BUTTONS.slider_step,'value',get(handles.BUTTONS.slider_step,'max'));
                        else
                            set(handles.BUTTONS.slider_step,'value',offset);
                        end
                        xlim = [offset-round(0.5*length_window) offset+round(0.5*length_window)];
                        xlim(1) = max([xlim(1) 1]);
                        xlim(2) = min([xlim(2) D.Nsamples]);

                        D.PSD.VIZU.xlim = xlim;

                        updateDisp(D)

                        set(handles.BUTTONS.focus_temp,'string',offset);
                        set(handles.BUTTONS.slider_step,'value',offset);

                    end

                end




                % Execute actions accessible from the event contextmenu : delete event
            case 'deleteEvent'

                D.trials.events(currentEvent) = [];

                D.PSD.tools.redo.select = D.trials.events;

                handles = rmfield(D.PSD.handles,'PLOT');
                D.PSD.handles = handles;
                updateDisp(D)

                if isempty(D.trials.events)
                    set(handles.SELECT.select_minus,'enable','off');
                    set(handles.SELECT.show_select,'enable','off');
                    set(handles.SELECT.save_select,'enable','off');
                    set(handles.TOOLS.cor_average,'enable','off');
                    set(handles.TOOLS.find_peaks,'enable','off');
                    set(handles.SELECT.select_nothing,'enable','off');
                    set(handles.SELECT.goto_select1,'enable','off');
                    set(handles.SELECT.goto_select2,'enable','off');
                    set(handles.TOOLS.spectrum_events,'enable','off');
                    set(handles.TOOLS.spectrum_comp,'enable','off');

                end


                %                 set(handles.EDIT.undo,'enable','on');
                %                 set(handles.EDIT.redo,'enable','off');

        end





        %% Selection callbacks
    case 'select'

        switch arg2


            %% Switch to another trial
            case 'switch'

                trN = get(gco,'value');
                if strcmp(D.PSD.VIZU.type,'scalpData')
                    handles = rmfield(D.PSD.handles,'PLOT');
                    D.PSD.handles = handles;
                end
                D.PSD.trials.current = trN;
                updateDisp(D)




                %% Add an event to current selection
            case 'add'

                D.PSD.tools.coreg      = 0;
                D.PSD.tools.undo.coreg = 0;

                [x,y]                             = ginput(2);
                x                               = round(x);
                x(1)                            = min([max([1 x(1)]) D.Nsamples]);
                x(2)                            = min([max([1 x(2)]) D.Nsamples]);
                x                               = sort(x(:)');
                Nevents = length(D.trials.events);
                D.trials.events(Nevents+1).time     = min(x)./D.Fsample;
                D.trials.events(Nevents+1).duration   = abs(diff(x))./D.Fsample;
                D.trials.events(Nevents+1).type       = '0';
                D.trials.events(Nevents+1).value      = 0;

                D.PSD.tools.redo.select        = D.trials.events;
                D.PSD.tools.redo.coreg         = 0;

                % Enable tools on selections
%                 set(handles.SELECT.select_minus,'enable','on');
%                 set(handles.SELECT.show_select,'enable','on');
%                 set(handles.SELECT.save_select,'enable','on');
%                 set(handles.TOOLS.cor_average,'enable','on');
%                 set(handles.TOOLS.find_peaks,'enable','on');
%                 set(handles.TOOLS.classify_peaks ,'enable','on');
%                 set(handles.SELECT.select_nothing,'enable','on');
%                 set(handles.SELECT.goto_select1,'enable','on');
%                 set(handles.SELECT.goto_select2,'enable','on');
                set(handles.BUTTONS.sb2,'enable','on');
                set(handles.BUTTONS.sb3,'enable','on');
%                 set(handles.TOOLS.spectrum_events,'enable','on');
                %     set(handles.TOOLS.spectrum_comp,'enable','on');
                %                 set(handles.EDIT.undo,'enable','on');
                %                 set(handles.EDIT.redo,'enable','off');

                % Update display
                handles = rmfield(D.PSD.handles,'PLOT');
                D.PSD.handles = handles;
                updateDisp(D)



                %% Remove last event in list
            case 'remove'

                D.PSD.tools.undo.coreg     = 0;
                D.PSD.tools.redo.select    = D.trials.events;
                D.PSD.tools.redo.coreg     = 0;

                D.PSD.tools.coreg          = 0;
                D.trials.events(:,end)        = [];

                % Disable tools
                if isempty(ud.select)
                    set(gcbo,'enable','off');
                    set(handles.SELECT.show_select,'enable','off');
                    set(handles.SELECT.save_select,'enable','off');
                    set(handles.TOOLS.cor_average,'enable','off');
                    set(handles.TOOLS.find_peaks,'enable','off');
                    set(handles.TOOLS.classify_peaks,'enable','off');
                    set(handles.SELECT.select_nothing,'enable','off');
                    set(handles.SELECT.goto_select1,'enable','off');
                    set(handles.SELECT.goto_select2,'enable','off');
                    set(handles.BUTTONS.sb2,'enable','off');
                    set(handles.BUTTONS.sb3,'enable','off');
                    set(handles.TOOLS.spectrum_events,'enable','off');
                    %         set(handles.TOOLS.spectrum_comp,'enable','off');

                end

                set(handles.EDIT.undo,'enable','on');
                set(handles.EDIT.redo,'enable','off');

                % Update display
                handles = rmfield(D.PSD.handles,'PLOT');
                D.PSD.handles = handles;
                updateDisp(D)


                %% scroll through data upto next event
            case 'goto'


                here                    = get(handles.BUTTONS.slider_step,'value');
%                 Nevents                 = length(D.trials.events);
                x                       = [D.trials.events.time]';
                x(:,2)                  = [D.trials.events.duration]';
                x(:,2)                  = sum(x,2);
                xm = mean(x,2).*D.Fsample;
                if arg3 == 0
                    ind = find(xm > here+1);
                else
                    ind = find(xm < here-1);
                end
                if ~isempty(ind)

                    if arg3 == 1
                        offset          = round(max(xm(ind)));
                    else
                        offset          = round(min(xm(ind)));
                    end

%                     ylim                = D.PSD.VIZU.ylim;
                    xlim0               = get(handles.axes,'xlim');

                    if ~isequal(xlim0,[1 D.Nsamples])

                        length_window   = round(xlim0(2)-xlim0(1));
                        if offset < round(0.5*length_window)
                            offset      = round(0.5*length_window);
                            set(handles.BUTTONS.slider_step,'value',1);
                        elseif offset > D.Nsamples-round(0.5*length_window)
                            offset      = D.Nsamples-round(0.5*length_window)-1;
                            set(handles.BUTTONS.slider_step,'value',get(handles.BUTTONS.slider_step,'max'));
                        else
                            set(handles.BUTTONS.slider_step,'value',offset);
                        end
                        xlim            = [offset-round(0.5*length_window) offset+round(0.5*length_window)];
                        xlim(1)         = max([xlim(1) 1]);
                        xlim(2)         = min([xlim(2) D.Nsamples]);

                        D.PSD.VIZU.xlim    = xlim;

%                         set(handles.BUTTONS.focus_temp,'string',offset);
                        set(handles.BUTTONS.slider_step,'value',offset);

                        updateDisp(D)

                    end

                end


                %% Build new data matrix from concatenated events
            case 'show'

                Nevents             = length(D.trials.events);
%                 offset              = 1;
                data_show           = [];
                x                   = [D.trials.events.time]';
                x(:,2)              = [D.trials.events.duration]';
                x(:,2)              = sum(x,2);
                x                   = floor(x.*D.Fsample) +1;
                for i = 1:Nevents
                    data_show       = [data_show,D.data.y(:,x(i,1):1:x(i,2))];
                end

                D2                  = D;
                D2.trials.events    = [];
                D2.data.y           = data_show;

                spm_eeg_review(D2)




                %% Save current selection
            case 'save'


                select_data.date                = date;
                select_data.select              = D.trials.events;
                Nselect = length(select_data);
                x                   = [D.trials.events.time]';
                x(:,2)              = [D.trials.events.duration]';
                x(:,2)              = sum(x,2);
                x                   = floor(x.*D.Fsample) +1;
                for i=1:Nselect
                    select_data.select(i).data  = D.data.y(:,x(i,1):1:x(i,2));
                end
                clear D handles


                uisave


                %% Load selection
            case 'load'


                button                              = questdlg(...
                    'This will erase all current selections. Are you sure?');

                if isequal(button,'Yes')
                    [filename,pathname]             = uigetfile(...
                        '.mat','Please choose selection file!');

                    if ~isequal(filename,0)

                        filename                    = fullfile(pathname,filename);
                        s                           = load(filename);
                        fn                          = fieldnames(s);

                        if ismember('select_data',fn)

                            select_data             = getfield(s,'select_data');
                            D.trials.events         = select_data.select;
                            %                 fprintf(1,'Checking selections...')
                            %                 % check selections format...
                            %                 Nevents                 = length(D.trials.events);
                            %                 for i = 1:Nevents
                            %                     ud.select(i).x    = sort(ud.select(i).x(:)');
                            %                     if ~isfield(ud.select(i),'eventType') | isempty(ud.select(i).eventType)
                            %                         ud.select(i).eventType = 1;
                            %                     end
                            %                 end
                            %                 fprintf(1,' OK.')
                            %                 fprintf(1,'\n')

                            set(handles.SELECT.select_minus,'enable','on');
                            set(handles.SELECT.show_select,'enable','on');
                            set(handles.SELECT.save_select,'enable','on');
                            set(handles.TOOLS.cor_average,'enable','on');
                            set(handles.TOOLS.find_peaks,'enable','on');
                            set(handles.TOOLS.classify_peaks,'enable','on');
                            set(handles.SELECT.select_nothing,'enable','on');
                            set(handles.SELECT.goto_select1,'enable','on');
                            set(handles.SELECT.goto_select2,'enable','on');
                            set(handles.BUTTONS.sb2,'enable','on');
                            set(handles.BUTTONS.sb3,'enable','on');
                            set(handles.TOOLS.spectrum_events,'enable','on');
                            %                 set(handles.TOOLS.spectrum_comp,'enable','on');

                            set(handles.EDIT.undo,'enable','on');
                            set(handles.EDIT.redo,'enable','off');

                            D.PSD.tools.redo.select    = D.trials.events;
                            D.PSD.tools.coreg          = 0;

                            handles = rmfield(D.PSD.handles,'PLOT');
                            D.PSD.handles = handles;
                            updateDisp(D)


                        else

                            h                       = msgbox(...
                                ['The file you have provided does not contain any PRESELECTDATA selection file.'],...
                                'Data import');
                            uiwait(h);

                        end

                    end
                end



                %% Select all data file
            case 'all'


                button                      = questdlg(...
                    'This will erase all current selections. Are you sure you want to replace your current selections by the whole data window?');

                if isequal(button,'Yes')

                    D.trials.events = [];
                    D.trials.events(1).time = 1./D.Fsample;
                    D.trials.events(1).duration = D.Nsamples.*D.Fsample;
                    D.trials.events(1).type = '0';
                    D.trials.events(1).value = 0;


                    set(handles.SELECT.select_minus,'enable','on');
                    set(handles.SELECT.show_select,'enable','on');
                    set(handles.SELECT.save_select,'enable','on');
                    set(handles.TOOLS.cor_average,'enable','on');
                    set(handles.TOOLS.find_peaks,'enable','on');
                    set(handles.SELECT.select_nothing,'enable','on');
                    set(handles.TOOLS.spectrum_events,'enable','on');
                    %     set(handles.TOOLS.spectrum_comp,'enable','on');
                    %         set(handles.SELECT.goto_select,'enable','on');

                    set(handles.EDIT.undo,'enable','on');
                    set(handles.EDIT.redo,'enable','off');

                    D.PSD.tools.redo.select    = D.trials.events;
                    D.PSD.tools.coreg          = 0;

                    handles = rmfield(D.PSD.handles,'PLOT');
                    D.PSD.handles = handles;
                    updateDisplay(D)

                end


                %% Erase all events in selection
            case 'nothing'

                button = questdlg('This will erase all current selections. Are you sure?');

                if isequal(button,'Yes')

                    D.trials.events = [];

                    set(handles.SELECT.select_minus,'enable','off');
                    set(handles.SELECT.show_select,'enable','off');
                    set(handles.SELECT.save_select,'enable','off');
                    set(handles.TOOLS.cor_average,'enable','off');
                    set(handles.TOOLS.find_peaks,'enable','off');
                    set(handles.TOOLS.classify_peaks,'enable','off');
                    set(handles.SELECT.select_nothing,'enable','off');
                    set(handles.SELECT.goto_select1,'enable','off');
                    set(handles.SELECT.goto_select2,'enable','off');
                    set(handles.BUTTONS.sb2,'enable','off');
                    set(handles.BUTTONS.sb3,'enable','off');
                    set(handles.TOOLS.spectrum_events,'enable','off');
                    %         set(handles.TOOLS.spectrum_comp,'enable','off');

                    D.PSD.tools.redo.select    = D.trials.events;
                    D.PSD.tools.redo.coreg     = D.PSD.tools.coreg;
                    D.PSD.tools.redo.coreg     = 0;

                    set(handles.EDIT.undo,'enable','on');
                    set(handles.EDIT.redo,'enable','off');

                    handles = rmfield(D.PSD.handles,'PLOT');
                    D.PSD.handles = handles;
                    updateDisplay(D)

                end

                %%
        end





        %% Tools callbacks
    case 'tools'


end



%% Main update display
function [] = updateDisp(D,flag)
% This function updates the display of the data and events.

dbstop if error
if ~exist('flag','var')
    flag = 0;
end
handles = D.PSD.handles;



% Create intermediary display variables : events
figure(handles.hfig)


switch D.PSD.VIZU.type

    case 'standardData'

        % Create intermediary display variables
        % Get [xlim,ylim] axes limits
%         if isequal(get(handles.VIZU.zoom2,'enable'),'on')
%             ylim                    = get(handles.axes(1),'ylim');
%         else
%             ylim                    = D.PSD.VIZU.ylim;
%         end
        ylim                        = D.PSD.VIZU.ylim;
        ylim0                       = D.PSD.VIZU.ylim0;
        xlim                        = sort(D.PSD.VIZU.xlim);
        nc = size(D.PSD.VIZU.montage.M,1);
%         inv = setdiff(1:nc,D.PSD.VIZU.visuSensors);
%         dx                          = abs(diff(xlim));
%         xs                          = max(xlim) - dx./10;
%         xs2                         = max(xlim) - dx./10.5;
%         dy                          = abs(diff(ylim));
%         ys                          = min(ylim) + dy./20;

        % Get data matrix and events to display
        if strcmp(D.PSD.type,'continuous') && ~isempty(D.trials.events)
            trN = 1;
            Nevents                 = length(D.trials.events);
%             col                     = colormap(lines);
%             col                     = col(1:7,:);
            x                       = [D.trials.events.time]';
            x(:,2)                  = [D.trials.events.duration]';
            x(:,2)                  = sum(x,2);
            x                       = x*D.Fsample;
            LookEvents              = find((x(:,1) <= xlim(2) & x(:,1) >= xlim(1))...
                | (x(:,2) <= xlim(2) & x(:,2) >= xlim(1))...
                | (x(:,1) <= xlim(2) & x(:,2) >= xlim(2)) );
            BlindEvents             = setdiff(1:Nevents,LookEvents);
        elseif strcmp(D.PSD.type,'epoched')
            trN = D.PSD.trials.current(1);
            Nevents = 0;
        elseif isempty(D.trials.events)
            trN = 1;
            Nevents   = 0;
        end
        
        v_data                  = full(D.PSD.VIZU.montage.M)*D.data.y(:,xlim(1):xlim(2),trN);
        v_data                  = D.PSD.VIZU.visu_scale*(v_data);

        % Create graphical objects if absent
        if ~isfield(handles,'PLOT')
            set(handles.axes,'xlim',xlim,'nextplot','add');
%             axes(handles.axes)
            % create uicontextmnu on channel time series
            % plot data on visualization window and add colour repairs on
            % window
            v_data = v_data +repmat(D.PSD.VIZU.offset,1,size(v_data,2));
%             handles.PLOT.p = plot(handles.axes,xlim(1):xlim(2),v_data');
%             hold on
%             handles.PLOT.p2 = plot(handles.axes,xlim(1),D.PSD.VIZU.offset,'s',...
%                 'markersize',2,'linewidth',4);
            col = colormap('lines');
            col = repmat(col(1:7,:),floor(nc./7)+1,1);
            for i=1:nc
                cmenu = uicontextmenu;
                uimenu(cmenu,'Label',['channel ',num2str(D.PSD.VIZU.visuSensors(i)),': ',D.PSD.VIZU.montage.clab{i}]);
                uimenu(cmenu,'Label',['type: ',D.channels(D.PSD.VIZU.visuSensors(i)).type]);
                uimenu(cmenu,'Label',['bad: ',num2str(D.channels(D.PSD.VIZU.visuSensors(i)).bad)],...
                    'callback',@switchBC,'userdata',i,...
                    'BusyAction','cancel',...
                    'Interruptible','off');
                status = D.channels(D.PSD.VIZU.visuSensors(i)).bad;
                if ~status
                    lineStyle = '-';
                else
                    lineStyle = ':';
                end
                handles.PLOT.p(i) = plot(handles.axes,xlim(1):xlim(2),v_data(i,:)',...
                    'uicontextmenu',cmenu,'lineStyle',lineStyle,...
                    'color',col(i,:));
                handles.PLOT.p2(i) = plot(handles.axes,xlim(1),D.PSD.VIZU.offset(i),'s',...
                    'markersize',2,'linewidth',4,...
                    'uicontextmenu',cmenu,...
                    'color',col(i,:));
%                 set(handles.PLOT.p(i),'uicontextmenu',cmenu,'lineStyle',lineStyle);
%                 set(handles.PLOT.p2(i),'uicontextmenu',cmenu);
            end
%             % set invisible channels
%             set(handles.PLOT.p(inv),'visible','off');
            % Add on patches/lines for visualization of events
            if Nevents > 0
                col                     = colormap(lines);
                col                     = col(1:7,:);
                values                  = [D.trials(trN).events.value];
                values                  = mod(values,7);
                values(values==0)       = 7;
                handles.PLOT.e          = zeros(Nevents,1);
                for i = 1:Nevents
                    if abs(diff(x(i,:))) >0 % create patch rectangle...
                        handles.PLOT.e(i)   = patch([x(i,1) x(i,1) x(i,2) x(i,2)],...
                            [ylim0(1) ylim0(2) ylim0(2) ylim0(1)],col(values(i),:),...
                            'parent',handles.axes);
                        set(handles.PLOT.e(i),'edgecolor','none','facealpha',0.30,...
                            'userdata',i,'ButtonDownFcn','set(gco,''selected'',''on'')');
                    else  % ... as well as left line marker (onset)
                        handles.PLOT.e(i)   = plot(handles.axes,[x(i,1) x(i,1)],...
                            [ylim0(1) ylim0(2)]);
                        set(handles.PLOT.e(i),'color',col(values(i),:),...
                            'userdata',i,'ButtonDownFcn','set(gco,''selected'',''on'')');
                    end
                    sc.currentEvent = i;
                    sc.eventType    = D.trials(trN).events(i).type;
                    sc.eventValue   = D.trials(trN).events(i).value;
                    sc.N_select     = Nevents;
                    psd_defineMenuEvent(handles.PLOT.e(i),sc);
                    if ismember(i,BlindEvents)
                        set(handles.PLOT.e(i),'visible','off')
                    end
                end
            end
            % Update axes limits and channel names
            hold off
            D.PSD.handles = handles;
            set(handles.axes,'xlim',xlim,'ylim',ylim,'ytick',D.PSD.VIZU.offset,...
                'yticklabel',D.PSD.VIZU.montage.clab,'fontsize',D.PSD.VIZU.FontSize);
            set(handles.hfig,'userdata',D);
        else
            v_data = v_data +repmat(D.PSD.VIZU.offset,1,size(v_data,2));
            % scroll through data
            for i=1:length(D.PSD.VIZU.visuSensors)
                set(handles.PLOT.p(i),...
                    'xdata',xlim(1):xlim(2),...
                    'ydata',v_data(i,:));
                set(handles.PLOT.p2(i),...
                    'xdata',xlim(1));
            end
%             % set visible channels
%             set(handles.PLOT.p(inv),'visible','off');
%             set(handles.PLOT.p(D.PSD.VIZU.visuSensors),'visible','on');
%             set(handles.PLOT.p2(inv),'visible','off');
%             set(handles.PLOT.p2(D.PSD.VIZU.visuSensors),'visible','on');
            % Add on patches for visualization of selected events
            if Nevents >0
                set(handles.PLOT.e(BlindEvents),'visible','off')
                set(handles.PLOT.e(LookEvents),'visible','on')
            end
            % Update axes limits and channel names
            set(handles.axes,'xlim',xlim)
            if ~flag
                set(handles.axes,'ylim',ylim,'ytick',D.PSD.VIZU.offset,...
                    'yticklabel',D.PSD.VIZU.montage.clab,'fontsize',D.PSD.VIZU.FontSize);
                set(handles.hfig,'userdata',D);
            end
        end
        % Update scale axes
        pos0 = get(handles.axes,'position');
        pos1 = get(handles.scale,'position');
        dt = (abs(diff(get(handles.axes,'xlim')))./D.Fsample).*(pos1(3)./pos0(3));
        dz = (abs(diff(get(handles.axes,'ylim')))).*(pos1(4)./pos0(4))./D.PSD.VIZU.visu_scale;
        set(handles.scale,'xticklabel',[num2str(dt.*1e3),' ms'],...
            'yticklabel',num2str(dz));


    case 'scalpData'
        
        if strcmp(D.transform.ID,'time')

            trN = D.PSD.trials.current;
            Ntrials = length(trN);
            v_data = [];
            for i=1:Ntrials
                v_datai                 = full(D.PSD.VIZU.montage.M)*D.data.y(:,:,trN(i));
                v_datai                 = D.PSD.VIZU.visu_scale*(v_datai);
                v_data(:,:,i)           = v_datai;
            end

            % Create graphical objects if absent
            if ~isfield(handles,'PLOT')

                miY = min(v_data(:));
                maY = max(v_data(:));

                for i=1:length(D.PSD.VIZU.visuSensors)
                    cmenu = uicontextmenu;
                    uimenu(cmenu,'Label',['channel ',num2str(D.PSD.VIZU.visuSensors(i)),': ',D.PSD.VIZU.montage.clab{i}]);
                    uimenu(cmenu,'Label',['type: ',D.channels(D.PSD.VIZU.visuSensors(i)).type]);
                    uimenu(cmenu,'Label',['bad: ',num2str(D.channels(D.PSD.VIZU.visuSensors(i)).bad)],...
                        'callback',@switchBC,'userdata',i,...
                        'BusyAction','cancel',...
                        'Interruptible','off');
                    status = D.channels(D.PSD.VIZU.visuSensors(i)).bad;
                    if ~status
                        color = [1 1 1];
                    else
                        color = 0.75*[1 1 1];
                    end

                    set(handles.fra(i),'uicontextmenu',cmenu);
                    set(handles.axes(i),'color',color,...
                        'ylim',[miY maY]./D.PSD.VIZU.visu_scale);
                    handles.PLOT.p(:,i) = plot(handles.axes(i),squeeze(v_data(i,:,:)),...
                        'uicontextmenu',cmenu,'userdata',i);

                end
                % Update axes limits and channel names
                D.PSD.handles = handles;

            else
                % scroll through data
                for i=1:length(D.PSD.VIZU.visuSensors)
                    for j=1:Ntrials
                        set(handles.PLOT.p(j,i),'ydata',v_data(i,:,j));
                    end
                end

            end
            % Update scale axes
            dz = (abs(diff(get(handles.axes(1),'ylim'))))./D.PSD.VIZU.visu_scale;
            set(handles.scale,'yticklabel',num2str(dz));
            set(handles.hfig,'userdata',D);
            axes(D.PSD.handles.scale)

        else %---- Time-frequency data !! ----%

            trN = D.PSD.trials.current;
            Ntrials = length(trN);

            miY = 0;
            maY = 0;
            
            for i=1:length(D.PSD.VIZU.visuSensors)

                datai = squeeze(D.data.y(i,:,:,trN(1)));
                miY = min([min(datai(:)),miY]);
                maY = max([max(datai(:)),maY]);
%                 if ~isfield(handles,'PLOT')
                    D.PSD.handles.PLOT.im(i) = imagesc(datai);
                    set(D.PSD.handles.PLOT.im(i),...
                        'parent',handles.axes(i),...
                        'userdata',i);
%                 else
%                     set(D.PSD.handles.PLOT.im(i),...
%                         'cdata',datai);
%                 end

            end
            for i=1:length(D.PSD.VIZU.visuSensors)
                caxis(handles.axes(i),[miY maY]);
                colormap('jet')
            end
            
            set(handles.hfig,'userdata',D);
            
            
        end
        
    case 'visuRecon'
        
        trN = D.PSD.trials.current(1);
        invN = D.PSD.invN;
        model = D.other.inv{invN}.inverse;

        J = zeros(model.Nd,size(model.T,1));
        J(model.Is,:) = model.J{trN}*model.T';
        
        time = (model.pst-D.PSD.VIZU.x0).^2;
        indTime = find(time==min(time));
        gridTime = model.pst(indTime);
        
        tex = J(:,indTime);

        set(D.PSD.handles.mesh,'facevertexcdata',tex)

        set(handles.hfig,'userdata',D);
        
        set(D.PSD.handles.BUTTONS.slider_step,'value',gridTime)
        set(D.PSD.handles.BUTTONS.focus_temp,'string',num2str(gridTime))

end




%% Switch 'bad channel' status
function [] = switchBC(varargin)
ind = get(gcbo,'userdata');
D = get(gcf,'userdata');
status = D.channels(D.PSD.VIZU.visuSensors(ind)).bad;
if status
    status = 0;
    lineStyle = '-';
    color = [1 1 1];
else
    status = 1;
    lineStyle = ':';
    color = 0.75*[1 1 1];
end
D.channels(D.PSD.VIZU.visuSensors(ind)).bad = status;
set(D.PSD.handles.hfig,'userdata',D);
cmenu = uicontextmenu;
uimenu(cmenu,'Label',['channel ',num2str(D.PSD.VIZU.visuSensors(ind)),': ',D.PSD.VIZU.montage.clab{ind}]);
uimenu(cmenu,'Label',['type: ',D.channels(D.PSD.VIZU.visuSensors(ind)).type]);
uimenu(cmenu,'Label',['bad: ',num2str(status)],...
    'callback',@switchBC,'userdata',ind,...
    'BusyAction','cancel',...
    'Interruptible','off');
switch D.PSD.VIZU.type
    case 'standardData'
        set(D.PSD.handles.PLOT.p(ind),'uicontextmenu',cmenu,...
            'lineStyle',lineStyle);
        set(D.PSD.handles.PLOT.p2(ind),'uicontextmenu',cmenu);
    case 'scalpData'
        set(D.PSD.handles.axes(ind),'Color',color);%,...
%             'uicontextmenu',cmenu);
        set(D.PSD.handles.fra(ind),'uicontextmenu',cmenu);
        set(D.PSD.handles.PLOT.p(:,ind),'uicontextmenu',cmenu);
        axes(D.PSD.handles.scale)
end



%% Define menu event
function [] = psd_defineMenuEvent(re,sc)
% This funcion defines the uicontextmenu associated to the selected events.
% All the actions which are accessible using the right mouse click on the
% selected events are a priori defined here.
dbstop if error
% Highlighting the selection
set(re,'buttondownfcn','spm_eeg_review_callbacks(''menuEvent'',''click'',0)');
cmenu = uicontextmenu;
set(re,'uicontextmenu',cmenu);
% Display basic info
info = ['--- EVENT #',num2str(sc.currentEvent),' /',...
    num2str(sc.N_select),' (type= ',sc.eventType,', value= ',num2str(sc.eventValue),') ---'];
uimenu(cmenu,'label',info,'enable','off');
% Properties editor
uimenu(cmenu,'separator','on','label','Edit event properties',...
    'callback','spm_eeg_review_callbacks(''menuEvent'',''EventProperties'',0)',...
    'BusyAction','cancel',...
    'Interruptible','off');
% Go to next event of the same type
hc = uimenu(cmenu,'label','Go to iso-type closest event');
uimenu(hc,'label','forward','callback','spm_eeg_review_callbacks(''menuEvent'',''goto'',1)',...
    'BusyAction','cancel',...
    'Interruptible','off');
uimenu(hc,'label','backward','callback','spm_eeg_review_callbacks(''menuEvent'',''goto'',0)',...
    'BusyAction','cancel',...
    'Interruptible','off');
% Delete action
uimenu(cmenu,'label','Delete event','callback','spm_eeg_review_callbacks(''menuEvent'',''deleteEvent'',0)',...
    'BusyAction','cancel',...
    'Interruptible','off');


%% Get info about source reconstruction
function str = getInfo4Inv(D,invN)
str{1} = ['Label: ',D.other.inv{invN}.comment{1}];
try
    str{2} = ['Date: ',D.other.inv{invN}.date(1,:),', ',D.other.inv{invN}.date(2,:)];
catch
    str{2} = ['Date: ',D.other.inv{invN}.date(1,:)];
end
str{3} = ['Modality: ',D.other.inv{invN}.modality];
if strcmp(D.other.inv{invN}.method,'Imaging')
    source = 'distributed';
else
    source = 'equivalent current dipoles';
end
str{4} = ['Source model: ',source,' (',D.other.inv{invN}.method,')'];
% str{5} = ['Nb of included dipoles: ',num2str(size(D.other.inv{invN}.inverse.J{1},1))];
str{5} = ['Nb of included dipoles: ',...
    num2str(length(D.other.inv{invN}.inverse.Is)),...
    ' / ',num2str(D.other.inv{invN}.inverse.Nd)];
str{6} = ['Inversion method: ',D.other.inv{invN}.inverse.type];
try
    str{7} = ['Time window of interest: ',...
        num2str(D.other.inv{invN}.inverse.woi(1)),...
        ' to ',num2str(D.other.inv{invN}.inverse.woi(2)),' ms'];
catch
    str{7} = ['Time window of interest: ',...
        num2str(D.other.inv{invN}.inverse.pst(1)),...
        ' to ',num2str(D.other.inv{invN}.inverse.pst(end)),' ms'];
end
try
    if D.other.inv{1}.inverse.Han
        han = 'yes';
    else
        han = 'no';
    end
    str{8} = ['Hanning: ',han];
catch
    str{8} = ['Hanning: ?'];
end
if isfield(D.other.inv{invN}.inverse,'lpf')
    str{9} = ['Band pass filter: ',num2str(D.other.inv{invN}.inverse.lpf),...
        ' to ',num2str(D.other.inv{invN}.inverse.hpf), 'Hz'];
else
    str{9} = ['Band pass filter: default'];
end
str{10} = ['Nb of temporal modes: ',...
    num2str(size(D.other.inv{invN}.inverse.T,2))];
str{11} = ['Variance accounted for: ',...
    num2str(D.other.inv{invN}.inverse.R2),' %'];
str{12} = ['Log model evidence (free energy): ',...
    num2str(D.other.inv{invN}.inverse.F)];


