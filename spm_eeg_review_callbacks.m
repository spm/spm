function [varargout] = spm_eeg_review_callbacks(varargin)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jean Daunizeau
% $Id: spm_eeg_review_callbacks.m 1979 2008-08-05 18:05:05Z jean $

try
    D = get(gcf,'userdata');
    handles = D.PSD.handles;
end

switch varargin{1}

    
    %% File I/O 
    case 'file'
        switch varargin{2}
            case 'save'
                spm('pointer','watch');
                drawnow
                D = rmfield(D,'PSD');
                D = meeg(D);
                D.save;
%                 save(D.fname,'D')
                spm('pointer','arrow');
        end
        
        
    case 'get'
        
        switch varargin{2}
            case 'VIZU'
                if strcmp(D.transform.ID,'time')
                    visuSensors             = varargin{3};
                    M                       = sparse(length(visuSensors),length(D.channels));
                    M(sub2ind(size(M),1:length(visuSensors),visuSensors(:)')) = 1;
                    decim                   = max([floor((D.Nsamples.*size(D.data.y,3))./1e4),1]);
                    data                    = D.data.y(visuSensors,1:decim:D.Nsamples,:);
                    sd                      = std(data(:));
                    offset                  = (0:1:length(visuSensors)-1)'*sd/2;
                    v_data                  = 0.25.*data +repmat(offset,[1 size(data,2) size(data,3)]);
                    ma                      = max(v_data(:))+sd;
                    mi                      = min(v_data(:))-sd;
                    ylim                    = [mi ma];
                    VIZU.visu_scale         = 0.25;
                    VIZU.FontSize           = 10;
                    VIZU.visuSensors        = visuSensors;
                    VIZU.visu_offset        = sd;
                    VIZU.offset             = offset;
                    VIZU.ylim               = ylim;
                    VIZU.ylim0              = ylim;
                    VIZU.figname            = 'main visualization window';
                    VIZU.montage.M          = M;
                    VIZU.montage.clab       = {D.channels(visuSensors).label};
                else
                    visuSensors             = varargin{3};
%                     VIZU.visuSensors     = visuSensors;
                    VIZU.montage.clab    = {D.channels(visuSensors).label};
                end
                varargout{1} = VIZU;
                return
            case 'commentInv'
                try
                    invN = varargin{3};
                catch
                    try
                        invN = D.PSD.invN;
                    catch
                        invN = 1;
                    end
                end
                str = getInfo4Inv(D,invN);
                varargout{1} = str;
                return
            case 'dataInfo'
                str = getInfo4Data(D);
                varargout{1} = str;
                return
            case 'uitable'
                D = getUItable(D);
                spm_eeg_review_switchDisplay(D);
            case 'prep'
                Finter = spm_figure('GetWin','Interactive');
                D = get(Finter, 'UserData');
                spm_eeg_review(D);
                spm_clf(Finter)
                
        end
    
    
    %% Visualization callbacks

    case 'visu'

        zoom off

        switch varargin{2}

            case 'main'

                switch varargin{3}
                    case 'eeg'
                        D.PSD.VIZU.modality = 'eeg';
                    case 'meg'
                        D.PSD.VIZU.modality = 'meg';
                    case 'other'
                        D.PSD.VIZU.modality = 'other';
                    case 'source'
                        D.PSD.VIZU.modality = 'source';
                    case 'info';
                        D.PSD.VIZU.modality = 'info';
                        try
                            D.PSD.VIZU.info = varargin{4};
                        catch
                            D.PSD.VIZU.info = 1;
                        end
                    case 'standard'
                        D.PSD.VIZU.type = 1;
                    case 'scalp'
                        D.PSD.VIZU.type = 2;
                end
                try
                    D.PSD.VIZU.xlim = get(handles.axes(1),'xlim');
                end
                [D] = spm_eeg_review_switchDisplay(D);
                try
                    updateDisp(D)
                end
                   

            case 'switch'
                
                spm('pointer','watch');
                drawnow
                mod = get(gcbo,'userdata');
                if ~isequal(mod,D.PSD.VIZU.type)
                    if mod == 1
                        spm_eeg_review_callbacks('visu','main','standard')
                    else
                        spm_eeg_review_callbacks('visu','main','scalp')
                    end
                end
                spm('pointer','arrow');
                
            case 'update'
                
                try
                    D = varargin{3};
                end
                updateDisp(D)
                
                % Scalp interpolation
            case 'scalp_interp'

                if ~isempty([D.channels(:).X_plot2D])
                    x = round(mean(get(handles.axes(1),'xlim')));
                    ylim = get(handles.axes(1),'ylim');
                    if D.PSD.VIZU.type==1
                        hl = line('parent',handles.axes,'xdata',[x;x],'ydata',[ylim(1);ylim(2)]);
                        in.hl = hl;
                    end
                    switch D.PSD.type
                        case 'continuous'
                            trN = 1;
                            in.gridTime = (1:D.Nsamples)./D.Fsample + D.timeOnset;
                            in.unit = 's';
                        case 'epoched'
                            trN = D.PSD.trials.current(1);
                            in.trN = trN;
                            in.gridTime = (1:D.Nsamples).*1e3./D.Fsample + D.timeOnset.*1e3;
                            in.unit = 'ms';
                    end
                    in.x = x;

                    in.handles = handles;
                    
                    switch D.PSD.VIZU.modality
                        case 'eeg'
                            I = D.PSD.EEG.I;
                            in.type = 'EEG';
                        case 'meg'
                            I = D.PSD.MEG.I;
                            in.type = 'MEG';
                        case 'other'
                            I = D.PSD.other.I;
                            in.type = 'other';
                    end
                    I = intersect(I,find(~[D.channels.bad]));
                    try
                        pos(:,1) = [D.channels(I).X_plot2D]';
                        pos(:,2) = [D.channels(I).Y_plot2D]';
                        labels = {D.channels(I).label};
                        y = D.data.y(I,:,trN);
                        in.min = min(y(:));
                        in.max = max(y(:));
                        in.ind = I;
                        in.type = 'EEG';
                        y = y(:,x);
                        spm_eeg_plotScalpData(y,pos,labels,in);
                    catch
                        msgbox('Get 2d positions for these channels!')
                    end
                else
                    msgbox('Get 2d positions for EEG/MEG channels!')
                end

            case 'inv'
                
%                 D.PSD.invN = varargin{3};
%                 if ~isequal(D.PSD.invN,varargin{3})
                    delete(findobj('tag','dipSpheres'))
                    str = getInfo4Inv(D,varargin{3});
                    isInv = get(D.PSD.handles.BMCcurrent,'userdata');
                    set(D.PSD.handles.infoText,'string',str);
                    set(D.PSD.handles.BMCcurrent,'XData',find(isInv==varargin{3}));
                    D.PSD.invN = varargin{3};
                    trN = D.PSD.trials.current(1);
                    model = D.other.inv{D.PSD.invN}.inverse;
                    J = model.J{trN}*model.T';
                    set(D.PSD.handles.axes,'CLim',[min(min(J)) max(max(J))]);
                    set(D.PSD.handles.mesh,...
                        'Vertices',D.other.inv{D.PSD.invN}.mesh.tess_mni.vert,...
                        'Faces',D.other.inv{D.PSD.invN}.mesh.tess_mni.face);
                    if isfield(D.other.inv{D.PSD.invN}.inverse,'dipfit') ||...
                            ~isequal(D.other.inv{D.PSD.invN}.inverse.xyz,zeros(1,3))
                        try 
                            xyz = D.other.inv{D.PSD.invN}.inverse.dipfit.Lpos;
                            radius = D.other.inv{D.PSD.invN}.inverse.dipfit.radius;
                        catch
                            xyz = D.other.inv{D.PSD.invN}.inverse.xyz';
                            radius = D.other.inv{D.PSD.invN}.inverse.rad(1);
                        end 
                        Np  = size(xyz,2);
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
%                 end


                % Contrast/intensity rescaling
            case 'iten_sc'

                switch D.PSD.VIZU.modality
                    case 'eeg'
                        D.PSD.EEG.VIZU.visu_scale = varargin{3}*D.PSD.EEG.VIZU.visu_scale;
                    case 'meg'
                        D.PSD.MEG.VIZU.visu_scale = varargin{3}*D.PSD.MEG.VIZU.visu_scale;
                    case 'other'
                        D.PSD.other.VIZU.visu_scale = varargin{3}*D.PSD.other.VIZU.visu_scale;
                end
                if D.PSD.VIZU.type==1
                    D.PSD.VIZU.xlim = get(D.PSD.handles.axes(1),'xlim');
                end
                updateDisp(D);


                % Resize plotted data window
            case 'time_w'

                % Get current plotted data window range and limits
                xlim = get(handles.axes(1),'xlim');
                %                 ylim = get(handles.axes(1),'ylim');
                length_window = max([varargin{3}*round(xlim(2)-xlim(1)),1]);
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
                if varargin{3} > 1
                    set(handles.BUTTONS.slider_step,'visible','on');
                end
                if ~isequal(xlim,[1 D.Nsamples])
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

                if varargin{3} > 1
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

                    case 1

                        if varargin{3}
                            zoom
                            %                     set(handles.VIZU.zoom2,'enable','on')
                        else % reset zoom and rebuild normal plotted data window
                            %                     set(handles.VIZU.zoom2,'enable','off')
                            updateDisp(D)
                        end

                    case 2
                        
                        switch D.PSD.VIZU.modality
                            case 'eeg'
                                VIZU = D.PSD.EEG.VIZU;
                            case 'meg'
                                VIZU = D.PSD.MEG.VIZU;
                            case 'other'
                                VIZU = D.PSD.other.VIZU;
                        end
                        
                        try
                            axes(D.PSD.handles.scale)
                        end
                        [x] = ginput(1);
                        indAxes = get(gco,'userdata');
                        if ~~indAxes
                            hf = figure;
                            chanLabel = D.channels(VIZU.visuSensors(indAxes)).label;
                            if D.channels(VIZU.visuSensors(indAxes)).bad
                                chanLabel = [chanLabel,' (BAD)'];
                            end
                            set(hf,'name',['channel ',chanLabel])
                            ha2 = axes('parent',hf,...
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
                                pst = (0:1./D.Fsample:(D.Nsamples-1)./D.Fsample) + D.timeOnset;
                                pst = pst*1e3;  % in msec
                                for i=1:Ntrials
                                    datai = get(hp(Ntrials-i+1),'ydata')./VIZU.visu_scale;
                                    plot(ha2,pst,datai,'color',col(i,:));
                                    leg{i} = D.PSD.trials.TrLabels{trN(i)};
                                end
                                legend(leg)
                                set(ha2,'xlim',[min(pst),max(pst)])
                                xlabel(ha2,'time (in ms after time onset)')
                                title(ha2,['channel ',chanLabel,...
                                    ' (',D.channels(VIZU.visuSensors(indAxes)).type,')'])

                            else % time-frequency data

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

                switch D.PSD.VIZU.modality
                    case 'eeg'
                        D.PSD.EEG.VIZU.visuSensors = D.PSD.EEG.I;
                    case 'meg'
                        D.PSD.MEG.VIZU.visuSensors = D.PSD.MEG.I;
                    case 'other'
                        D.PSD.other.VIZU.visuSensors = D.PSD.other.I;
                end
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
                if ~strcmp(D.PSD.VIZU.modality,'source')
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
                        %                         set(handles.BUTTONS.focus
                        %                         _temp,'string',round(mean(xlim)));
                        D.PSD.VIZU.xlim = xlim;
                        updateDisp(D,1)
                    end
                else
                    D.PSD.source.VIZU.x0 = offset;
                   
                    updateDisp(D,1)
                end

                %% Scroll page by page (button)
            case 'goOne'

                % Get current plotted data window range and limits
                xlim0 = get(handles.axes(1),'xlim');
                xm = mean(xlim0);
                length_window = abs(diff(xlim0));
                if varargin{3} == 0
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

                switch D.PSD.VIZU.modality
                    case 'eeg'
                        D.PSD.EEG.VIZU.montage.M = -D.PSD.EEG.VIZU.montage.M;
                    case 'meg'
                        D.PSD.MEG.VIZU.montage.M = -D.PSD.MEG.VIZU.montage.M;
                    case 'other'
                        D.D.PSD.other.VIZU.montage.M = -D.PSD.other.VIZU.montage.M;
                end
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


        switch varargin{2}

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
                        eventValue = str2double(answer{2});
                        D.trials.events(currentEvent).time = str2double(answer{3});
                        D.trials.events(currentEvent).duration = str2double(answer{4});
                        D.trials.events(currentEvent).type = eventType;
                        D.trials.events(currentEvent).value = eventValue;
                    end

                    D.PSD.tools.redo.select = D.trials.events;

                    try
                        delete(D.PSD.handles.PLOT.p)
                    end
                    try
                        delete(D.PSD.handles.PLOT.p2)
                    end
                    try
                        delete(D.PSD.handles.PLOT.e)
                    end
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
                if varargin{3} == 0
                    ind = find(xm < here);
                else
                    ind = find(xm > here);
                end

                if ~isempty(ind)

                    if varargin{3} == 0
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

                try
                    delete(D.PSD.handles.PLOT.p)
                end
                try
                    delete(D.PSD.handles.PLOT.p2)
                end
                try
                    delete(D.PSD.handles.PLOT.e)
                end
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

        switch varargin{2}


            %% Switch to another trial
            case 'switch'

                trN = get(gco,'value');
                if D.PSD.VIZU.type == 2
                    handles = rmfield(D.PSD.handles,'PLOT');
                    D.PSD.handles = handles;
                end
                D.PSD.trials.current = trN;
                status = [prod([D.trials(trN).bad])];
                try
                    if status
                        str = ['declare as not bad'];
                    else
                        str = ['declare as bad'];
                    end
                    set(D.PSD.handles.BUTTONS.badEvent,'string',str);
                end
                updateDisp(D)


            case 'bad'
                
                trN = D.PSD.trials.current;
                str = get(D.PSD.handles.BUTTONS.badEvent,'string');
                str1 = 'not bad';
                str2 = 'bad';
                if strcmp(str,['declare as ',str2])
                    bad = 1;
                    lab = [' (',str2,')'];
                    str = ['declare as ',str1];
                else
                    bad = 0;
                    lab = [' (',str1,')'];
                    str = ['declare as ',str2];
                end
                nt = length(trN);
                for i=1:nt
                    D.trials(trN(i)).bad = bad;
                    D.PSD.trials.TrLabels{trN(i)} = ['Trial ',num2str(trN(i)),...
                        ': ',D.trials(trN(i)).label,lab];
                end
                set(D.PSD.handles.BUTTONS.pop1,'string',D.PSD.trials.TrLabels);
                set(D.PSD.handles.BUTTONS.badEvent,'string',str)
                set(D.PSD.handles.hfig,'userdata',D)
                

                %% Add an event to current selection
            case 'add'

                D.PSD.tools.coreg = 0;
                D.PSD.tools.undo.coreg = 0;

                [x] = ginput(2);
                x = round(x);
                x(1) = min([max([1 x(1)]) D.Nsamples]);
                x(2) = min([max([1 x(2)]) D.Nsamples]);
                x = sort(x(:)');
                Nevents = length(D.trials.events);
                D.trials.events(Nevents+1).time = min(x)./D.Fsample;
                D.trials.events(Nevents+1).duration = abs(diff(x))./D.Fsample;
                D.trials.events(Nevents+1).type = '0';
                D.trials.events(Nevents+1).value = 0;

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
                try
                    delete(D.PSD.handles.PLOT.p)
                end
                try
                    delete(D.PSD.handles.PLOT.p2)
                end
                try
                    delete(D.PSD.handles.PLOT.e)
                end
                handles = rmfield(D.PSD.handles,'PLOT');
                D.PSD.handles = handles;
                updateDisp(D)



                %% Remove last event in list
            case 'remove'

                D.PSD.tools.undo.coreg     = 0;
                D.PSD.tools.redo.select    = D.trials.events;
                D.PSD.tools.redo.coreg     = 0;
                D.PSD.tools.coreg          = 0;
                D.trials.events(:,end)     = [];
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
                try
                    delete(D.PSD.handles.PLOT.p)
                end
                try
                    delete(D.PSD.handles.PLOT.p2)
                end
                try
                    delete(D.PSD.handles.PLOT.e)
                end
                handles = rmfield(D.PSD.handles,'PLOT');
                D.PSD.handles = handles;
                updateDisp(D)


                %% scroll through data upto next event
            case 'goto'


                here                    = get(handles.BUTTONS.slider_step,'value');
                x                       = [D.trials.events.time]';
                x(:,2)                  = [D.trials.events.duration]';
                x(:,2)                  = sum(x,2);
                xm = mean(x,2).*D.Fsample;
                if varargin{3} == 0
                    ind = find(xm > here+1);
                else
                    ind = find(xm < here-1);
                end
                if ~isempty(ind)
                    if varargin{3} == 1
                        offset          = round(max(xm(ind)));
                    else
                        offset          = round(min(xm(ind)));
                    end
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

                            try
                                delete(D.PSD.handles.PLOT.p)
                            end
                            try
                                delete(D.PSD.handles.PLOT.p2)
                            end
                            try
                                delete(D.PSD.handles.PLOT.e)
                            end
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

                    try
                        delete(D.PSD.handles.PLOT.p)
                    end
                    try
                        delete(D.PSD.handles.PLOT.p2)
                    end
                    try
                        delete(D.PSD.handles.PLOT.e)
                    end
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

                    try
                        delete(D.PSD.handles.PLOT.p)
                    end
                    try
                        delete(D.PSD.handles.PLOT.p2)
                    end
                    try
                        delete(D.PSD.handles.PLOT.e)
                    end
                    handles = rmfield(D.PSD.handles,'PLOT');
                    D.PSD.handles = handles;
                    updateDisplay(D)

                end

                %%
        end

        %% Edit callbacks (from spm_eeg_prep_ui)
    case 'edit'
        
        switch varargin{2}
            
            case 'prep'
                spm_eeg_prep_ui;
                Finter = spm_figure('GetWin','Interactive');
                D = rmfield(D,'PSD');
                if isempty(D.other)
                    D.other = 1;
                end
                D.other.PSD = 1;
                D = meeg(D);
                set(Finter, 'UserData', D);
                hc = get(Finter,'children');
                delete(hc(end));    % get rid of 'file' uimenu...
                %... and add an 'OK' button:
                uicontrol(Finter,...
                    'style','pushbutton','string','OK',...
                    'callback','spm_eeg_review_callbacks(''get'',''prep'')',...
                    'tooltipstring','Update data informations in ''SPM Graphics'' window',...
                    'BusyAction','cancel',...
                    'Interruptible','off',...
                    'Tag','EEGprepUI');
                set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'File'), 'Enable', 'on');
                IsEEG = 'off';
                IsMEG = 'off';
                HasSensors = 'off';
                HasSensorsEEG = 'off';
                HasSensorsMEG = 'off';
                Dloaded = 'on';
                if ~isempty(strmatch('EEG', D.chantype, 'exact'))
                    IsEEG = 'on';
                end
                if ~isempty(strmatch('MEG', D.chantype, 'exact'));
                    IsMEG = 'on';
                end
                if ~isempty(D.sensors('EEG')) || ~isempty(D.sensors('MEG'))
                    HasSensors = 'on';
                end
                if ~isempty(D.sensors('EEG'))
                    HasSensorsEEG = 'on';
                end
                if  ~isempty(D.sensors('MEG'))
                    HasSensorsMEG = 'on';
                end
                IsTemplate = 'off';
                IsSelected = 'off';
                IsMoved = 'off';
                set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Channel types'), 'Enable', Dloaded);
                set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Sensors'), 'Enable', Dloaded);
                set(findobj(Finter,'Tag','EEGprepUI', 'Label', '2D projection'), 'Enable', Dloaded);
                set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Load EEG sensors'), 'Enable', IsEEG);
                set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Load MEG Fiducials/Headshape'), 'Enable', HasSensorsMEG);
                set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Headshape'), 'Enable', HasSensors);
                set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Coregister (EEG)'), 'Enable', HasSensorsEEG);
                set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Coregister (MEG)'), 'Enable', HasSensorsMEG);
                set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Edit existing EEG'), 'Enable', IsEEG);
                set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Edit existing MEG'), 'Enable', IsMEG);
                set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Project 3D (EEG)'), 'Enable', HasSensorsEEG);
                set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Project 3D (MEG)'), 'Enable', HasSensorsMEG);
                set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Delete sensor'), 'Enable', IsSelected);
                set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Undo move'), 'Enable', IsMoved);
                set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Apply'), 'Enable', IsTemplate);
                set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Clear'), 'Enable', IsTemplate);
                delete(setdiff(findobj(Finter), [Finter; findobj(Finter,'Tag','EEGprepUI')]));
                figure(Finter);

            
            case 'sensorPos'
                
                spm_eeg_prep_ui('LoadEEGSensCB')
                spm_eeg_prep_ui('HeadshapeCB')
            case 'coregSpos'
                spm_eeg_prep_ui('CoregisterCB')
            case 'volSens'
                spm_eeg_prep_ui('PrepVolSensCB')
            case '2Dpos'
                spm_eeg_prep_ui('EditExistingCoor2DCB')
            case 'project3Dpos'
                spm_eeg_prep_ui('Project3DCB')
                
                
                
        end

        %% Tools callbacks
    case 'tools'


end



%% Main update display
function [] = updateDisp(D,flag)
% This function updates the display of the data and events.

if ~exist('flag','var')
    flag = 0;
end
handles = D.PSD.handles;



% Create intermediary display variables : events
figure(handles.hfig)

if ~strcmp(D.PSD.VIZU.modality,'source')

    switch D.PSD.VIZU.modality
        case 'eeg'
            VIZU = D.PSD.EEG.VIZU;
        case 'meg'
            VIZU = D.PSD.MEG.VIZU;
        case 'other'
            VIZU = D.PSD.other.VIZU;
    end

    switch D.PSD.VIZU.type

        case 1
            
            % Create intermediary display variables
            ylim                        = VIZU.ylim;
            ylim0                       = VIZU.ylim0;
            xlim                        = sort(D.PSD.VIZU.xlim);
            nc = size(VIZU.montage.M,1);
            % Get data matrix and events to display
            if strcmp(D.PSD.type,'continuous') && ~isempty(D.trials.events)
                trN = 1;
                Nevents                 = length(D.trials.events);
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
            v_data                  = full(VIZU.montage.M)*D.data.y(:,xlim(1):xlim(2),trN);
            v_data                  = VIZU.visu_scale*(v_data);

            % Create graphical objects if absent
            if ~isfield(handles,'PLOT')
                set(handles.axes,'xlim',xlim,'nextplot','add');
                % create uicontextmnu on channel time series
                % plot data on visualization window and add colour repairs on
                % window
                v_data = v_data +repmat(VIZU.offset,1,size(v_data,2));
                col = colormap('lines');
                col = repmat(col(1:7,:),floor(nc./7)+1,1);
                for i=1:nc
                    cmenu = uicontextmenu;
                    uimenu(cmenu,'Label',['channel ',num2str(VIZU.visuSensors(i)),': ',VIZU.montage.clab{i}]);
                    uimenu(cmenu,'Label',['type: ',D.channels(VIZU.visuSensors(i)).type]);
                    uimenu(cmenu,'Label',['bad: ',num2str(D.channels(VIZU.visuSensors(i)).bad)],...
                        'callback',@switchBC,'userdata',i,...
                        'BusyAction','cancel',...
                        'Interruptible','off');
                    status = D.channels(VIZU.visuSensors(i)).bad;
                    if ~status
                        lineStyle = '-';
                    else
                        lineStyle = ':';
                    end
                    handles.PLOT.p(i) = plot(handles.axes,xlim(1):xlim(2),v_data(i,:)',...
                        'uicontextmenu',cmenu,'lineStyle',lineStyle,...
                        'color',col(i,:),'tag','plotEEG');
                    handles.PLOT.p2(i) = plot(handles.axes,xlim(1),VIZU.offset(i),'s',...
                        'markersize',2,'linewidth',4,...
                        'uicontextmenu',cmenu,...
                        'color',col(i,:),'tag','plotEEG');
                end
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
                                'userdata',i,'ButtonDownFcn','set(gco,''selected'',''on'')',...
                                'tag','plotEEG');
                        else  % ... as well as left line marker (onset)
                            handles.PLOT.e(i)   = plot(handles.axes,[x(i,1) x(i,1)],...
                                [ylim0(1) ylim0(2)]);
                            set(handles.PLOT.e(i),'color',col(values(i),:),...
                                'userdata',i,'ButtonDownFcn','set(gco,''selected'',''on'')',...
                                'tag','plotEEG');
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
                set(handles.axes,'xlim',xlim,'ylim',ylim,'ytick',VIZU.offset,...
                    'yticklabel',VIZU.montage.clab,'fontsize',VIZU.FontSize);
                set(handles.hfig,'userdata',D);
            else
                v_data = v_data +repmat(VIZU.offset,1,size(v_data,2));
                % scroll through data
                for i=1:length(VIZU.visuSensors)
                    set(handles.PLOT.p(i),...
                        'xdata',xlim(1):xlim(2),...
                        'ydata',v_data(i,:));
                    set(handles.PLOT.p2(i),...
                        'xdata',xlim(1));
                end
                % Add on patches for visualization of selected events
                if Nevents >0
                    set(handles.PLOT.e(BlindEvents),'visible','off')
                    set(handles.PLOT.e(LookEvents),'visible','on')
                end
                % Update axes limits and channel names
                set(handles.axes,'xlim',xlim)
                if ~flag
                    set(handles.axes,'ylim',ylim,'ytick',VIZU.offset,...
                        'yticklabel',VIZU.montage.clab,'fontsize',VIZU.FontSize);
                    set(handles.hfig,'userdata',D);
                end
            end
            % Update scale axes
            pos0 = get(handles.axes,'position');
            pos1 = get(handles.scale,'position');
            dt = (abs(diff(get(handles.axes,'xlim')))./D.Fsample).*(pos1(3)./pos0(3));
            dz = (abs(diff(get(handles.axes,'ylim')))).*(pos1(4)./pos0(4))./VIZU.visu_scale;
            set(handles.scale,'xticklabel',[num2str(dt.*1e3),' ms'],...
                'yticklabel',num2str(dz));


        case 2

            if strcmp(D.transform.ID,'time')
                
                trN = D.PSD.trials.current;
                Ntrials = length(trN);
                v_data = zeros(size(VIZU.montage.M,1),...
                    size(D.data.y,2),Ntrials);
%                 v_data = [];
                for i=1:Ntrials
                    v_datai                 = full(VIZU.montage.M)*D.data.y(:,:,trN(i));
                    v_datai                 = VIZU.visu_scale*(v_datai);
                    v_data(:,:,i)           = v_datai;
                end
                % Create graphical objects if absent
                if ~isfield(handles,'PLOT')
                    miY = min(v_data(:));
                    maY = max(v_data(:));
                    for i=1:length(VIZU.visuSensors)
                        cmenu = uicontextmenu;
                        uimenu(cmenu,'Label',['channel ',num2str(VIZU.visuSensors(i)),': ',VIZU.montage.clab{i}]);
                        uimenu(cmenu,'Label',['type: ',D.channels(VIZU.visuSensors(i)).type]);
                        uimenu(cmenu,'Label',['bad: ',num2str(D.channels(VIZU.visuSensors(i)).bad)],...
                            'callback',@switchBC,'userdata',i,...
                            'BusyAction','cancel',...
                            'Interruptible','off');
                        status = D.channels(VIZU.visuSensors(i)).bad;
                        if ~status
                            color = [1 1 1];
                        else
                            color = 0.75*[1 1 1];
                        end

                        set(handles.fra(i),'uicontextmenu',cmenu);
                        set(handles.axes(i),'color',color,...
                            'ylim',[miY maY]./VIZU.visu_scale);
                        handles.PLOT.p(:,i) = plot(handles.axes(i),squeeze(v_data(i,:,:)),...
                            'uicontextmenu',cmenu,'userdata',i,'tag','plotEEG');
                    end
                    % Update axes limits and channel names
                    D.PSD.handles = handles;
                else
                    % scroll through data
                    for i=1:length(VIZU.visuSensors)
                        for j=1:Ntrials
                            set(handles.PLOT.p(j,i),'ydata',v_data(i,:,j));
                        end
                    end
                end
                % Update scale axes
                dz = (abs(diff(get(handles.axes(1),'ylim'))))./VIZU.visu_scale;
                set(handles.scale,'yticklabel',num2str(dz));
                set(handles.hfig,'userdata',D);
                axes(D.PSD.handles.scale)

            else %---- Time-frequency data !! ----%

                trN = D.PSD.trials.current;
                miY = 0;
                maY = 0;
                for i=1:length(VIZU.visuSensors)
                    datai = squeeze(D.data.y(i,:,:,trN(1)));
                    miY = min([min(datai(:)),miY]);
                    maY = max([max(datai(:)),maY]);
                    D.PSD.handles.PLOT.im(i) = imagesc(datai);
                    set(D.PSD.handles.PLOT.im(i),'tag',plotEEG');
                    set(D.PSD.handles.PLOT.im(i),...
                        'parent',handles.axes(i),...
                        'userdata',i);
                end
                for i=1:length(VIZU.visuSensors)
                    caxis(handles.axes(i),[miY maY]);
                    colormap('jet')
                end
                set(handles.hfig,'userdata',D);

            end
    end
    
    
else  % source space

    trN = D.PSD.trials.current(1);
    invN = D.PSD.invN;
    model = D.other.inv{invN}.inverse;

    J = zeros(model.Nd,size(model.T,1));
    J(model.Is,:) = model.J{trN}*model.T';
    
    time = (model.pst-D.PSD.source.VIZU.x0).^2;
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
switch D.PSD.VIZU.modality
    case 'eeg'
        I = D.PSD.EEG.I;
        VIZU = D.PSD.EEG.VIZU;
    case 'meg'
        I = D.PSD.MEG.I;
        VIZU = D.PSD.MEG.VIZU;
    case 'other'
        I = D.PSD.other.I;
        VIZU = D.PSD.other.VIZU;
end
status = D.channels(I(ind)).bad;
if status
    status = 0;
    lineStyle = '-';
    color = [1 1 1];
else
    status = 1;
    lineStyle = ':';
    color = 0.75*[1 1 1];
end
D.channels(I(ind)).bad = status;
set(D.PSD.handles.hfig,'userdata',D);
cmenu = uicontextmenu;
uimenu(cmenu,'Label',['channel ',num2str(I(ind)),': ',VIZU.montage.clab{ind}]);
uimenu(cmenu,'Label',['type: ',D.channels(I(ind)).type]);
uimenu(cmenu,'Label',['bad: ',num2str(status)],...
    'callback',@switchBC,'userdata',ind,...
    'BusyAction','cancel',...
    'Interruptible','off');
switch D.PSD.VIZU.type
    case 1
        set(D.PSD.handles.PLOT.p(ind),'uicontextmenu',cmenu,...
            'lineStyle',lineStyle);
        set(D.PSD.handles.PLOT.p2(ind),'uicontextmenu',cmenu);
    case 2
        set(D.PSD.handles.axes(ind),'Color',color);
        set(D.PSD.handles.fra(ind),'uicontextmenu',cmenu);
        set(D.PSD.handles.PLOT.p(:,ind),'uicontextmenu',cmenu);
        axes(D.PSD.handles.scale)
end



%% Define menu event
function [] = psd_defineMenuEvent(re,sc)
% This funcion defines the uicontextmenu associated to the selected events.
% All the actions which are accessible using the right mouse click on the
% selected events are a priori defined here.

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
str{5} = ['Nb of included dipoles: ',...
    num2str(length(D.other.inv{invN}.inverse.Is)),...
    ' / ',num2str(D.other.inv{invN}.inverse.Nd)];
str{6} = ['Inversion method: ',D.other.inv{invN}.inverse.type];
try
    str{7} = ['Time window: ',...
        num2str(floor(D.other.inv{invN}.inverse.woi(1))),...
        ' to ',num2str(floor(D.other.inv{invN}.inverse.woi(2))),' ms'];
catch
    str{7} = ['Time window: ',...
        num2str(floor(D.other.inv{invN}.inverse.pst(1))),...
        ' to ',num2str(floor(D.other.inv{invN}.inverse.pst(end))),' ms'];
end
try
    if D.other.inv{invN}.inverse.Han
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


%% Get data info
function str = getInfo4Data(D)
str{1} = ['File name: ',D.path,filesep,D.fname];
str{2} = ['Type: ',D.type];
if ~strcmp(D.transform.ID,'time')
    str{2} = [str{2},' (time-frequency data)'];
end
str{3} = ['Number of samples: ',num2str(D.Nsamples)];
str{4} = ['Sampling frequency: ',num2str(D.Fsample)];
nb = length(find([D.channels.bad]));
str{5} = ['Number of channels: ',num2str(length(D.channels)),' (',num2str(nb),' bad channels)'];
nb = length(find([D.trials.bad]));
if strcmp(D.type,'continuous')
    str{6} = ['Number of events: ',num2str(length(D.trials(1).events))];
else
    str{6} = ['Number of trials: ',num2str(length(D.trials)),' (',num2str(nb),' bad trials)'];
end
try
    str{7} = ['Time onset: ',num2str(D.timeOnset)];
end


%% extracting data from spm_uitable java object
function [D] = getUItable(D)
spm('pointer','watch');
drawnow
ht = D.PSD.handles.infoUItable;
cn = get(ht,'columnNames');
table = get(ht,'data');
if length(cn) == 5  % channel info
    nc = length(D.channels);
    for i=1:nc
        if ~isempty(table(i,1))
            D.channels(i).label = table(i,1);
        end
        if ~isempty(table(i,2))
            switch lower(table(i,2))
                case 'eeg'
                    D.channels(i).type = 'EEG';
                case 'meg'
                    D.channels(i).type = 'MEG';
                case 'lfp'
                    D.channels(i).type = 'LFP';
                case 'veog'
                    D.channels(i).type = 'VEOG';
                case 'heog'
                    D.channels(i).type = 'HEOG';
                case 'other'
                    D.channels(i).type = 'Other';
                otherwise
                    D.channels(i).type = 'Other';
            end
        end
        if ~isempty(table(i,3))
            switch lower(table(i,3))
                case 'yes'
                    D.channels(i).bad = 1;
                otherwise
                    D.channels(i).bad = 0;
            end
        end
        if ~isempty(table(i,5))
            D.channels(i).units = table(i,5);
        end
    end
    % Find indices of channel types (these might have been changed)
    D.PSD.EEG.I  = find(strcmp('EEG',{D.channels.type}));
    D.PSD.MEG.I  = find(strcmp('MEG',{D.channels.type}));
    D.PSD.other.I = setdiff(1:nc,[D.PSD.EEG.I(:);D.PSD.MEG.I(:)]);
    if ~isempty(D.PSD.EEG.I)
        [out] = spm_eeg_review_callbacks('get','VIZU',D.PSD.EEG.I);
        D.PSD.EEG.VIZU = out;
    end
    if ~isempty(D.PSD.MEG.I)
        [out] = spm_eeg_review_callbacks('get','VIZU',D.PSD.MEG.I);
        D.PSD.MEG.VIZU = out;
    end
    if ~isempty(D.PSD.other.I)
        [out] = spm_eeg_review_callbacks('get','VIZU',D.PSD.other.I);
        D.PSD.other.VIZU = out;
    end
elseif length(cn) == 7
    if strcmp(D.type,'continuous')
        ne = length(D.trials(1).events);
        D.trials = rmfield(D.trials,'events');
        j = 0;
        for i=1:ne
            if isempty(table(i,1))&&...
                    isempty(table(i,2))&&...
                    isempty(table(i,3))&&...
                    isempty(table(i,4))&&...
                    isempty(table(i,5))&&...
                    isempty(table(i,6))&&...
                    isempty(table(i,7))
                % Row (ie event) has been cleared/deleted
            else
                j = j+1;
                if ~isempty(table(i,2))
                    D.trials(1).events(j).type = table(i,2);
                end
                if ~isempty(table(i,3))
                    D.trials(1).events(j).value = str2double(table(i,3));
                end
                if ~isempty(table(i,4))
                    D.trials(1).events(j).duration = str2double(table(i,4));
                end
                if ~isempty(table(i,5))
                    D.trials(1).events(j).time = str2double(table(i,5));
                end
            end
        end
    else
        nt = length(D.trials);
        for i=1:nt
            if ~isempty(table(i,1))
                D.trials(i).label = table(i,1);
            end
            if ~isempty(table(i,2))
                D.trials(i).events.type = table(i,2);
            end
            if ~isempty(table(i,3))
                D.trials(i).events.value = str2double(table(i,3));
            end
            if ~isempty(table(i,6))
                switch lower(table(i,6))
                    case 'yes'
                        D.trials(i).bad = 1;
                    otherwise
                        D.trials(i).bad = 0;
                end
            end
            if D.trials(i).bad
                str = ' (bad)';
            else
                str = ' (not bad)';
            end
            D.PSD.trials.TrLabels{i} = ['Trial ',num2str(i),': ',D.trials(i).label,str];
        end
    end

elseif length(cn) == 3
    nt = length(D.trials);
    for i=1:nt
        if ~isempty(table(i,1))
            D.trials(i).label = table(i,1);
        end
        if ~isempty(table(i,3))
            switch lower(table(i,3))
                case 'yes'
                    D.channels(i).bad = 1;
                otherwise
                    D.channels(i).bad = 0;
            end
        end
        if D.trials(i).bad
            str = ' (bad)';
        else
            str = ' (not bad)';
        end
        D.PSD.trials.TrLabels{i} = ['Trial ',num2str(i),': ',D.trials(i).label,str];
    end
    
elseif length(cn) == 12     % source reconstructions
    if isfield(D.other,'inv') && ~isempty(D.other.inv)
        inv = D.other.inv;
        Ninv = length(inv);
        D.other = rmfield(D.other,'inv');
        j = 0;
        for i=1:Ninv
            if isempty(table(i,1))&&...
                    isempty(table(i,2))&&...
                    isempty(table(i,3))&&...
                    isempty(table(i,4))&&...
                    isempty(table(i,5))&&...
                    isempty(table(i,6))&&...
                    isempty(table(i,7))&&...
                    isempty(table(i,8))&&...
                    isempty(table(i,9))&&...
                    isempty(table(i,10))&&...
                    isempty(table(i,11))&&...
                    isempty(table(i,12))
                % Row (ie source reconstruction) has been cleared/deleted
            else
                j = j+1;
                D.other.inv{j} = inv{i};
            end
        end
    end
end
set(D.PSD.handles.hfig,'userdata',D)
spm_eeg_review_callbacks('visu','main','info',D.PSD.VIZU.info)
drawnow
spm('pointer','arrow');

