function [varargout] = spm_eeg_review_callbacks(varargin)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jean Daunizeau
% $Id: spm_eeg_review_callbacks.m 2284 2008-10-01 15:54:05Z jean $

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
                spm('pointer','arrow');
        end


    case 'get'

        switch varargin{2}
            case 'VIZU'
                if strcmp(D.transform.ID,'time')
                    visuSensors             = varargin{3};
                    M                       = sparse(length(visuSensors),length(D.channels));
                    M(sub2ind(size(M),1:length(visuSensors),visuSensors(:)')) = 1;
                    if isequal(D.type,'continuous')
                        decim                   = max([floor((D.Nsamples.*size(D.data.y,3))./2e2),1]);
                    else
                        decim = 1;
                    end
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
                    VIZU.y2                 = permute(sum(data.^2,1),[2 3 1]);
                    VIZU.sci                = size(VIZU.y2,1)./D.Nsamples;
                else
                    visuSensors             = varargin{3};
                    VIZU.visuSensors        = visuSensors;
                    VIZU.montage.clab       = {D.channels(visuSensors).label};
                end
                varargout{1} = VIZU;
                return
            case 'commentInv'
                invN = varargin{3};
                str = getInfo4Inv(D,invN);
                varargout{1} = str;
                return
            case 'dataInfo'
                str = getInfo4Data(D);
                varargout{1} = str;
                return
            case 'history'
                table = getHistory(D);
                varargout{1} = table;
                return
            case 'uitable'
                D = getUItable(D);
                spm_eeg_review_switchDisplay(D);
            case 'prep'
                Finter = spm_figure('GetWin','Interactive');
                D = struct(get(Finter, 'UserData'));
                other = rmfield(D.other,'PSD');
                D.other = other;
                spm_eeg_review(D);
                spm_clf(Finter)

        end


        %% Visualization callbacks

    case 'visu'

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
                        end
                    case 'standard'
                        D.PSD.VIZU.type = 1;
                    case 'scalp'
                        D.PSD.VIZU.type = 2;
                end
                try,D.PSD.VIZU.xlim = get(handles.axes(1),'xlim');end
                [D] = spm_eeg_review_switchDisplay(D);
                try,updateDisp(D);end


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

                try,D = varargin{3};end
                updateDisp(D)

                % Scalp interpolation
            case 'scalp_interp'

                if ~isempty([D.channels(:).X_plot2D])
                    x = round(mean(get(handles.axes(1),'xlim')));
                    ylim = get(handles.axes(1),'ylim');
                    if D.PSD.VIZU.type==1
                        hl = line('parent',handles.axes,'xdata',[x;x],...
                            'ydata',[ylim(1);ylim(2)]);
                        in.hl = hl;
                    end
                    switch D.PSD.type
                        case 'continuous'
                            trN = 1;
                        case 'epoched'
                            trN = D.PSD.trials.current(1);
                            in.trN = trN;
                    end
                    in.gridTime = (1:D.Nsamples).*1e3./D.Fsample + D.timeOnset.*1e3;
                    in.unit = 'ms';
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


            case 'sensorPos'

                % get 3D positions
                try     % EEG
                    pos3d = [D.sensors.eeg.pnt];
                    m.vertices = D.other.inv{1}.mesh.tess_mni.vert;
                    m.faces = D.other.inv{1}.mesh.tess_mni.face;
                    options.figname = 'EEG sensors';
                    [out] = spm_eeg_render(m,options);
                    lo = load(fullfile(spm('Dir'),'EEGtemplates','wmeshTemplate_scalp'));
                    m2.vertices = lo.vert;
                    m2.faces = lo.face;
                    options.hfig = out.handles.fi;
                    options.ParentAxes = gca;
                    [out] = spm_eeg_render(m2,options);
                    figure(out.handles.fi)
                    hold on
                    plot3(pos3d(:,1),pos3d(:,2),pos3d(:,3),'.');
                    %                     text(pos3d(:,1),pos3d(:,2),pos3d(:,3),D.PSD.EEG.VIZU.montage.clab);
                    text(pos3d(:,1),pos3d(:,2),pos3d(:,3),D.sensors.eeg.label);
                    axis tight
                end
                try     % MEG
                    pos3d = [D.sensors.meg.pnt];
                    m.vertices = D.other.inv{1}.mesh.tess_mni.vert;
                    m.faces = D.other.inv{1}.mesh.tess_mni.face;
                    options.figname = 'MEG sensors';
                    [out] = spm_eeg_render(m,options);
                    lo = load(fullfile(spm('Dir'),'EEGtemplates','wmeshTemplate_scalp'));
                    m2.vertices = lo.vert;
                    m2.faces = lo.face;
                    options.hfig = out.handles.fi;
                    options.ParentAxes = gca;
                    [out] = spm_eeg_render(m2,options);
                    figure(out.handles.fi)
                    hold on
                    plot3(pos3d(:,1),pos3d(:,2),pos3d(:,3),'.');
                    %                     text(pos3d(:,1),pos3d(:,2),pos3d(:,3),D.PSD.MEG.VIZU.montage.clab);
                    text(pos3d(:,1),pos3d(:,2),pos3d(:,3),D.sensors.meg.label);
                    axis tight
                end
                try     % other
                    pos3d = [D.sensors.other.pnt];
                    m.vertices = D.other.inv{1}.mesh.tess_mni.vert;
                    m.faces = D.other.inv{1}.mesh.tess_mni.face;
                    options.figname = 'other sensors';
                    [out] = spm_eeg_render(m,options);
                    lo = load(fullfile(spm('Dir'),'EEGtemplates','wmeshTemplate_scalp'));
                    m2.vertices = lo.vert;
                    m2.faces = lo.face;
                    options.hfig = out.handles.fi;
                    options.ParentAxes = gca;
                    [out] = spm_eeg_render(m2,options);
                    figure(out.handles.fi)
                    hold on
                    plot3(pos3d(:,1),pos3d(:,2),pos3d(:,3),'.');
                    %                     text(pos3d(:,1),pos3d(:,2),pos3d(:,3),D.PSD.other.VIZU.montage.clab);
                    ht = text(pos3d(:,1),pos3d(:,2),pos3d(:,3),D.sensors.other.label);
                    axis tight
                end

            case 'inv'

                cla(D.PSD.handles.axes2,'reset')
                D.PSD.source.VIZU.current = varargin{3};
                updateDisp(D);


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
                updateDisp(D,1);


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
                    set(handles.BUTTONS.vb3,'enable','off')
                elseif length_window < 20
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
                    set(handles.BUTTONS.vb4,'enable','on');
                else
                    set(handles.BUTTONS.vb3,'enable','on');
                end

                updateDisp(D,1)


                %% Data navigation using the slider
            case 'slider_t'

                offset = get(gco,'value');
                if ~strcmp(D.PSD.VIZU.modality,'source')
                    offset = round(offset);
                    % Get current plotted data window range and limits
                    xlim0 = get(handles.axes(1),'xlim');
                    % The IF statement ensures acceptable range
                    if ~isequal(xlim0,[1 D.Nsamples])
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
                        updateDisp(D)
                    end

                else
                    updateDisp(D)

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

                % Zoom (box in)
            case 'zoom'

                switch D.PSD.VIZU.type

                    case 1

                        if ~isempty(D.PSD.handles.zoomh)
                            switch get(D.PSD.handles.zoomh,'enable')
                                case 'on'
                                    set(D.PSD.handles.zoomh,'enable','off')
                                case 'off'
                                    set(D.PSD.handles.zoomh,'enable','on')
                            end
                        else
                            if get(D.PSD.handles.BUTTONS.vb5,'value')
                                zoom on;
                            else
                                zoom off;
                            end
                            %set(D.PSD.handles.BUTTONS.vb5,'value',~val);
                        end

                    case 2

                        set(D.PSD.handles.BUTTONS.vb5,'value',1)
                        switch D.PSD.VIZU.modality
                            case 'eeg'
                                VIZU = D.PSD.EEG.VIZU;
                            case 'meg'
                                VIZU = D.PSD.MEG.VIZU;
                            case 'other'
                                VIZU = D.PSD.other.VIZU;
                        end
                        try,axes(D.PSD.handles.scale);end
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
                                'nextplot','add',...
                                'XGrid','on','YGrid','on');
                            trN = D.PSD.trials.current(:);
                            Ntrials = length(trN);

                            if strcmp(D.transform.ID,'time')

                                leg = cell(Ntrials,1);
                                col = colormap('lines');
                                col = repmat(col(1:7,:),floor(Ntrials./7)+1,1);
                                hp = get(handles.axes(indAxes),'children');
                                pst = (0:1/D.Fsample:(D.Nsamples-1)/D.Fsample) + D.timeOnset;
                                pst = pst*1e3;  % in msec
                                for i=1:Ntrials
                                    datai = get(hp(Ntrials-i+1),'ydata')./VIZU.visu_scale;
                                    plot(ha2,pst,datai,'color',col(i,:));
                                    leg{i} = D.PSD.trials.TrLabels{trN(i)};
                                end
                                legend(leg)
                                set(ha2,'xlim',[min(pst),max(pst)],...
                                    'ylim',get(D.PSD.handles.axes(indAxes),'ylim'))
                                xlabel(ha2,'time (in ms after time onset)')
                                title(ha2,['channel ',chanLabel,...
                                    ' (',D.channels(VIZU.visuSensors(indAxes)).type,')'])

                            else % time-frequency data

                                datai = squeeze(D.data.y(indAxes,:,:,trN(1)));
                                hp2 = image(datai,'CDataMapping','scaled');
                                set(hp2,'parent',ha2);
                                colormap('jet')
                                colorbar
                                pst = (0:1/D.Fsample:(D.Nsamples-1)/D.Fsample) + D.timeOnset;
                                pst = pst*1e3;  % in msec
                                set(ha2,'xtick',1:10:length(pst),'xticklabel',pst(1:10:length(pst)),...
                                    'xlim',[1 length(pst)],...
                                    'ytick',1:length(D.transform.frequencies),...
                                    'yticklabel',D.transform.frequencies);
                                xlabel(ha2,'time (in ms after time onset)')
                                ylabel(ha2,'frequency (in Hz)')
                                title(ha2,['channel ',chanLabel,...
                                    ' (',D.channels(VIZU.visuSensors(indAxes)).type,')'])

                            end

                            axes(ha2)
                        end
                        set(D.PSD.handles.BUTTONS.vb5,'value',0)
                end


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
                stc = cell(4,1);
                default = cell(4,1);
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

                    try,delete(D.PSD.handles.PLOT.p);end
                    try,delete(D.PSD.handles.PLOT.p2);end
                    try,delete(D.PSD.handles.PLOT.e);end
                    handles = rmfield(D.PSD.handles,'PLOT');
                    D.PSD.handles = handles;
                    updateDisp(D)

                end


                % Execute actions accessible from the event contextmenu : go to next/previous event
            case 'goto'


                here = mean(x(currentEvent,:));
                values = [D.trials.events.value];
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
                        set(handles.BUTTONS.slider_step,'value',offset);
                    end
                end



                % Execute actions accessible from the event contextmenu : delete event
            case 'deleteEvent'

                D.trials.events(currentEvent) = [];
                try,delete(D.PSD.handles.PLOT.p);end
                try,delete(D.PSD.handles.PLOT.p2);end
                try,delete(D.PSD.handles.PLOT.e);end
                handles = rmfield(D.PSD.handles,'PLOT');
                D.PSD.handles = handles;
                updateDisp(D)

        end





        %% Selection callbacks
    case 'select'

        switch varargin{2}


            %% Switch to another trial
            case 'switch'
                trN = get(gco,'value');
                if ~strcmp(D.PSD.VIZU.modality,'source') && D.PSD.VIZU.type == 2
                    handles = rmfield(D.PSD.handles,'PLOT');
                    D.PSD.handles = handles;
                else
                    try,cla(D.PSD.handles.axes2,'reset');end
                end
                D.PSD.trials.current = trN;
                status = ~any(~[D.trials(trN).bad]);
                try
                    if status
                        str = ['declare as not bad'];
                    else
                        str = ['declare as bad'];
                    end
                    ud = get(D.PSD.handles.BUTTONS.badEvent,'userdata');
                    set(D.PSD.handles.BUTTONS.badEvent,...
                        'tooltipstring',str,...
                        'cdata',ud.img{2-status},'userdata',ud)
                    switch D.PSD.VIZU.modality
                        case 'eeg'
                            VIZU = D.PSD.EEG.VIZU;
                        case 'meg'
                            VIZU = D.PSD.MEG.VIZU;
                        case 'other'
                            VIZU = D.PSD.other.VIZU;
                    end
                    set(handles.GV.y2,'ydata',VIZU.y2(:,trN));
                    set(handles.GV.axes,'ylim',[min(VIZU.y2(:,trN))...
                        max(VIZU.y2(:,trN))]);
                    set(handles.GV.p,'ydata',[min(VIZU.y2(:,trN)) ...
                        max(VIZU.y2(:,trN)) max(VIZU.y2(:,trN)) min(VIZU.y2(:,trN))]);
                end
                updateDisp(D)

            case 'bad'
                trN = D.PSD.trials.current;
                ud = get(D.PSD.handles.BUTTONS.badEvent,'userdata');
                str1 = 'not bad';
                str2 = 'bad';
                if ud.val
                    bad = 0;
                    lab = [' (',str1,')'];
                    str = ['declare as ',str2];
                else
                    bad = 1;
                    lab = [' (',str2,')'];
                    str = ['declare as ',str1];
                end
                ud.val = bad;
                nt = length(trN);
                for i=1:nt
                    D.trials(trN(i)).bad = bad;
                    D.PSD.trials.TrLabels{trN(i)} = ['Trial ',num2str(trN(i)),...
                        ': ',D.trials(trN(i)).label,lab];
                end
                set(D.PSD.handles.BUTTONS.pop1,'string',D.PSD.trials.TrLabels);
                set(D.PSD.handles.BUTTONS.badEvent,...
                    'tooltipstring',str,...
                    'cdata',ud.img{2-bad},'userdata',ud)
                set(D.PSD.handles.hfig,'userdata',D)
                try
                    uicontrol(D.PSD.handles.BUTTONS.pop1)
                end


                %                 str = get(D.PSD.handles.BUTTONS.badEvent,'string');
                %                 str1 = 'not bad';
                %                 str2 = 'bad';
                %                 if strcmp(str,['declare as ',str2])
                %                     bad = 1;
                %                     lab = [' (',str2,')'];
                %                     str = ['declare as ',str1];
                %                 else
                %                     bad = 0;
                %                     lab = [' (',str1,')'];
                %                     str = ['declare as ',str2];
                %                 end
                %                 nt = length(trN);
                %                 for i=1:nt
                %                     D.trials(trN(i)).bad = bad;
                %                     D.PSD.trials.TrLabels{trN(i)} = ['Trial ',num2str(trN(i)),...
                %                         ': ',D.trials(trN(i)).label,lab];
                %                 end
                %                 set(D.PSD.handles.BUTTONS.pop1,'string',D.PSD.trials.TrLabels);
                %                 set(D.PSD.handles.BUTTONS.badEvent,'string',str)
                %                 set(D.PSD.handles.hfig,'userdata',D)
                %                 try
                %                     uicontrol(D.PSD.handles.BUTTONS.pop1)
                %                 end

                %% Add an event to current selection
            case 'add'
                [x,tmp] = ginput(2);
                x = round(x);
                x(1) = min([max([1 x(1)]) D.Nsamples]);
                x(2) = min([max([1 x(2)]) D.Nsamples]);
                x = sort(x(:)');
                Nevents = length(D.trials.events);
                D.trials.events(Nevents+1).time = min(x)./D.Fsample;
                D.trials.events(Nevents+1).duration = abs(diff(x))./D.Fsample;
                D.trials.events(Nevents+1).type = '0';
                D.trials.events(Nevents+1).value = 0;
                % Enable tools on selections
                set(handles.BUTTONS.sb2,'enable','on');
                set(handles.BUTTONS.sb3,'enable','on');
                % Update display
                try,delete(D.PSD.handles.PLOT.p),end
                try,delete(D.PSD.handles.PLOT.p2),end
                try,delete(D.PSD.handles.PLOT.e),end
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
                        set(handles.BUTTONS.slider_step,'value',offset);
                        updateDisp(D)
                    end
                end

        end

        %% Edit callbacks (from spm_eeg_prep_ui)
    case 'edit'

        switch varargin{2}

            case 'prep'
                try,rotate3d off;end
                spm_eeg_prep_ui;
                Finter = spm_figure('GetWin','Interactive');
                D = rmfield(D,'PSD');
                if isempty(D.other)
                    D.other = struct([]);
                end
                D.other(1).PSD = 1;
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

                spm_eeg_prep_ui('update_menu')
                delete(setdiff(findobj(Finter), [Finter; findobj(Finter,'Tag','EEGprepUI')]));
                figure(Finter);

        end


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
            xlim                        = D.PSD.VIZU.xlim;
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
                set(handles.axes,'xlim',xlim,'ylim',ylim,'ytick',VIZU.offset,...
                    'yticklabel',VIZU.montage.clab,'fontsize',VIZU.FontSize);
                % Update scale axes
                pos0 = get(handles.axes,'position');
                pos1 = get(handles.scale,'position');
                dt = (abs(diff(get(handles.axes,'xlim')))./D.Fsample).*(pos1(3)./pos0(3));
                dz = (abs(diff(get(handles.axes,'ylim')))).*(pos1(4)./pos0(4))./VIZU.visu_scale;
                set(handles.scale,'xticklabel',[num2str(dt.*1e3),' ms'],...
                    'yticklabel',num2str(dz));
                % Display global data squared power
                handles.GV.y2 = plot(handles.GV.axes,VIZU.y2(:,trN));
                xlim2 = VIZU.sci.*[xlim(1) xlim(2)];
                if abs(diff(xlim2)) <= 1
                    xlim2(2) = xlim2(1) +1;
                end
                set(handles.GV.axes,'xlim',[1 size(VIZU.y2,1)],...
                    'ylim',[min(VIZU.y2(:,trN)) max(VIZU.y2(:,trN))],...
                    'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],...
                    'tag','plotEEG','box','on');
                handles.GV.p   = patch([xlim2(1) xlim2(1) xlim2(2) xlim2(2)],...
                    [min(VIZU.y2(:,trN)) max(VIZU.y2(:,trN)) max(VIZU.y2(:,trN)) min(VIZU.y2(:,trN))],...
                    [0.5 0.5 0.5],'parent',handles.GV.axes);
                set(handles.GV.p,'edgecolor','none','facealpha',0.8,...
                    'tag','plotEEG');
                D.PSD.handles = handles;
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
                % And update scale axes
                set(handles.axes,'xlim',xlim)
                if flag
                    set(handles.axes,'ylim',ylim,'ytick',VIZU.offset,...
                        'yticklabel',VIZU.montage.clab,'fontsize',VIZU.FontSize);
                    set(handles.hfig,'userdata',D);
                    pos0 = get(handles.axes,'position');
                    pos1 = get(handles.scale,'position');
                    dt = (abs(diff(get(handles.axes,'xlim')))./D.Fsample).*(pos1(3)./pos0(3));
                    dz = (abs(diff(get(handles.axes,'ylim')))).*(pos1(4)./pos0(4))./VIZU.visu_scale;
                    set(handles.scale,'xticklabel',[num2str(dt.*1e3),' ms'],...
                        'yticklabel',num2str(dz));
                end
                % Update position of global power patch
                xlim2 = VIZU.sci.*[xlim(1) xlim(2)];
                if abs(diff(xlim2)) <= 1
                    xlim2(2) = xlim2(1) +1;
                end
                set(handles.GV.p,'xdata',[xlim2(1) xlim2(1) xlim2(2) xlim2(2)])
            end


        case 2

            if strcmp(D.transform.ID,'time')

                trN = D.PSD.trials.current;
                Ntrials = length(trN);
                v_data = zeros(size(VIZU.montage.M,1),...
                    size(D.data.y,2),Ntrials);
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
                    cmenu = uicontextmenu;
                    uimenu(cmenu,'Label',['channel ',num2str(VIZU.visuSensors(i)),': ',VIZU.montage.clab{i}]);
                    uimenu(cmenu,'Label',['type: ',D.channels(VIZU.visuSensors(i)).type]);
                    %                     uimenu(cmenu,'Label',['bad: ',num2str(D.channels(VIZU.visuSensors(i)).bad)],...
                    %                         'callback',@switchBC,'userdata',i,...
                    %                         'BusyAction','cancel',...
                    %                         'Interruptible','off');
                    status = D.channels(VIZU.visuSensors(i)).bad;
                    if ~status
                        color = [1 1 1];
                    else
                        color = 0.75*[1 1 1];
                    end
                    datai = squeeze(D.data.y(i,:,:,trN(1)));
                    miY = min([min(datai(:)),miY]);
                    maY = max([max(datai(:)),maY]);
                    D.PSD.handles.PLOT.im(i) = image(datai,'CDataMapping','scaled');
                    set(D.PSD.handles.PLOT.im(i),...
                        'tag','plotEEG',...
                        'parent',handles.axes(i),...
                        'userdata',i,...
                        'hittest','off');
                    set(handles.fra(i),'uicontextmenu',cmenu);
                end
                colormap(jet)
                % This for normalized colorbars:
                %                 for i=1:length(VIZU.visuSensors)
                %                     caxis(handles.axes(i),[miY maY]);
                %                     colormap('jet')
                %                 end
                set(handles.hfig,'userdata',D);

            end
    end


else  % source space

    % get model/trial info
    VIZU = D.PSD.source.VIZU;
    invN = VIZU.isInv(D.PSD.source.VIZU.current);
    model = D.other.inv{invN}.inverse;
    t0 = get(D.PSD.handles.BUTTONS.slider_step,'value');
    tmp = (model.pst-t0).^2;
    indTime = find(tmp==min(tmp));
    gridTime = model.pst(indTime);

    try % simple time scroll
        % update time line
        set(VIZU.lineTime,'xdata',[gridTime;gridTime]);
        % update mesh's texture
        tex = VIZU.J(:,indTime);
        set(D.PSD.handles.mesh,'facevertexcdata',tex)
        set(D.PSD.handles.BUTTONS.slider_step,'value',gridTime)

    catch % VIZU.lineTime deleted -> switch to another source recon
        % get the inverse model info
        str = getInfo4Inv(D,invN);
        set(D.PSD.handles.infoText,'string',str);
        try, set(D.PSD.handles.BMCcurrent,'XData',invN); end;
        % get model/trial time series
        trN = D.PSD.trials.current(1);
        D.PSD.source.VIZU.J = zeros(model.Nd,size(model.T,1));
        D.PSD.source.VIZU.J(model.Is,:) = model.J{trN}*model.T';
        D.PSD.source.VIZU.miJ = min(min(D.PSD.source.VIZU.J));
        D.PSD.source.VIZU.maJ = max(max(D.PSD.source.VIZU.J));
        % modify mesh/texture and add spheres...
        tex = D.PSD.source.VIZU.J(:,indTime);
        set(D.PSD.handles.axes,'CLim',...
            [D.PSD.source.VIZU.miJ D.PSD.source.VIZU.maJ]);
        set(D.PSD.handles.mesh,...
            'Vertices',D.other.inv{invN}.mesh.tess_mni.vert,...
            'Faces',D.other.inv{invN}.mesh.tess_mni.face,...
            'facevertexcdata',tex);
        try; delete(D.PSD.handles.dipSpheres);end
        if isfield(D.other.inv{invN}.inverse,'dipfit') ||...
                ~isequal(D.other.inv{invN}.inverse.xyz,zeros(1,3))
            try
                xyz = D.other.inv{invN}.inverse.dipfit.Lpos;
                radius = D.other.inv{invN}.inverse.dipfit.radius;
            catch
                xyz = D.other.inv{invN}.inverse.xyz';
                radius = D.other.inv{invN}.inverse.rad(1);
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
        % modify time series plot itself
        switch D.PSD.source.VIZU.timeCourses
            case 1
                D.PSD.source.VIZU.plotTC = plot(D.PSD.handles.axes2,...
                    model.pst,D.PSD.source.VIZU.J');
            otherwise
                % this is meant to be extended for displaying something
                % else than just J (e.g. J^2, etc...)
        end
        grid(D.PSD.handles.axes2,'on')
        box(D.PSD.handles.axes2,'on')
        xlabel(D.PSD.handles.axes2,'peri-stimulus time (ms)')
        ylabel(D.PSD.handles.axes2,'sources intensity')
        % add time line repair
        set(D.PSD.handles.axes2,...
            'ylim',[D.PSD.source.VIZU.miJ,D.PSD.source.VIZU.maJ],...
            'xlim',[D.PSD.source.VIZU.pst(1),D.PSD.source.VIZU.pst(end)],...
            'nextplot','add');
        D.PSD.source.VIZU.lineTime = line('parent',D.PSD.handles.axes2,...
            'xdata',[gridTime;gridTime],...
            'ydata',[D.PSD.source.VIZU.miJ,D.PSD.source.VIZU.maJ]);
        set(D.PSD.handles.axes2,'nextplot','replace',...
            'tag','plotEEG');
        % change time slider value if out of bounds
        set(D.PSD.handles.BUTTONS.slider_step,'value',gridTime)
        % update data structure
        set(handles.hfig,'userdata',D);

    end


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
try
    str{5} = ['Nb of included dipoles: ',...
        num2str(length(D.other.inv{invN}.inverse.Is)),...
        ' / ',num2str(D.other.inv{invN}.inverse.Nd)];
catch
    str{5} = 'Nb of included dipoles: undefined';
end
try
    str{6} = ['Inversion method: ',D.other.inv{invN}.inverse.type];
catch
    str{6} = 'Inversion method: undefined';
end
try
    try
        str{7} = ['Time window: ',...
            num2str(floor(D.other.inv{invN}.inverse.woi(1))),...
            ' to ',num2str(floor(D.other.inv{invN}.inverse.woi(2))),' ms'];
    catch
        str{7} = ['Time window: ',...
            num2str(floor(D.other.inv{invN}.inverse.pst(1))),...
            ' to ',num2str(floor(D.other.inv{invN}.inverse.pst(end))),' ms'];
    end
catch
    str{7} = 'Time window: undefined';
end
try
    if D.other.inv{invN}.inverse.Han
        han = 'yes';
    else
        han = 'no';
    end
    str{8} = ['Hanning: ',han];
catch
    str{8} = ['Hanning: undefined'];
end
try
    if isfield(D.other.inv{invN}.inverse,'lpf')
        str{9} = ['Band pass filter: ',num2str(D.other.inv{invN}.inverse.lpf),...
            ' to ',num2str(D.other.inv{invN}.inverse.hpf), 'Hz'];
    else
        str{9} = ['Band pass filter: default'];
    end
catch
    str{9} = 'Band pass filter: undefined';
end
try
    str{10} = ['Nb of temporal modes: ',...
        num2str(size(D.other.inv{invN}.inverse.T,2))];
catch
    str{10} = 'Nb of temporal modes: undefined';
end
try
    str{11} = ['Variance accounted for: ',...
        num2str(D.other.inv{invN}.inverse.R2),' %'];
catch
    str{11} = 'Variance accounted for: undefined';
end
try
    str{12} = ['Log model evidence (free energy): ',...
        num2str(D.other.inv{invN}.inverse.F)];
catch
    str{12} = 'Log model evidence (free energy): undefined';
end


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
try,str{7} = ['Time onset: ',num2str(D.timeOnset)];end


%% Get history info
function table = getHistory(D)
try
    history = D.history;
    nf = length(history);
    table = cell(nf,3);
    for i=1:nf
        table{i,1} = history(i).fun;
        switch history(i).fun
            case 'spm_eeg_convert'
                table{i,2} = D.history(i).args.dataset;
                [path,fn] = fileparts(table{i,2});
                if i<nf
                    table{i,3} = fullfile(path,D.history(i).args.outfile);
                else
                    table{i,3} = '[this file]';
                end
            case 'spm_eeg_prep'
                table{i,1} = [table{i,1},' (',D.history(i).args.task,')'];
                try
                    [path,fn] = fileparts(table{i-1,3});
                catch
                    path = [];
                end
                table{i,2} = fullfile(path,D.history(i).args.D);
                if i<nf
                    table{i,3} = fullfile(path,D.history(i).args.D);
                else
                    table{i,3} = '[this file]';
                end
            otherwise
                table{i,2} = D.history(i).args.D;
                if i<nf
                    try
                        table{i,3} = D.history(i+1).args.D;
                    catch
                        table{i,3} = '?';
                    end
                else
                    table{i,3} = '[this file]';
                end
        end
    end
catch
    history = [];
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
            ne = length(D.trials(i).events);
            if ne<2
                if ~isempty(table(i,2))
                    D.trials(i).events.type = table(i,2);
                end
                if ~isempty(table(i,3))
                    D.trials(i).events.value = str2double(table(i,3));
                end
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
        D.PSD.trials.TrLabels{i} = ['Trial ',num2str(i),' (average of ',...
            num2str(D.trials(i).repl),' events): ',D.trials(i).label];
    end

elseif length(cn) == 12     % source reconstructions

    if ~~D.PSD.source.VIZU.current
        isInv = D.PSD.source.VIZU.isInv;
        inv = D.other.inv;
        Ninv = length(inv);
        D.other = rmfield(D.other,'inv');
        oV = D.PSD.source.VIZU;
        D.PSD.source = rmfield(D.PSD.source,'VIZU');
        pst = [];
        j = 0;  % counts the total number of final inverse solutions in D
        k = 0;  % counts the number of original 'imaging' inv sol
        l = 0;  % counts the number of final 'imaging' inv sol
        for i=1:Ninv
            if ~ismember(i,isInv)   % not 'imaging' inverse solutions
                j = j+1;
                D.other.inv{j} = inv{i};
            else                    % 'imaging' inverse solutions
                k = k+1;
                if isempty(table(k,1))&&...
                        isempty(table(k,2))&&...
                        isempty(table(k,3))&&...
                        isempty(table(k,4))&&...
                        isempty(table(k,5))&&...
                        isempty(table(k,6))&&...
                        isempty(table(k,7))&&...
                        isempty(table(k,8))&&...
                        isempty(table(k,9))&&...
                        isempty(table(k,10))&&...
                        isempty(table(k,11))&&...
                        isempty(table(k,12))
                    % Row (ie source reconstruction) has been cleared/deleted
                    % => erase inverse solution from D struct
                else
                    j = j+1;
                    l = l+1;
                    pst = [pst;inv{isInv(k)}.inverse.pst(:)];
                    D.other.inv{j} = inv{isInv(k)};
                    D.other.inv{j}.comment{1} = table(k,1);
                    D.PSD.source.VIZU.isInv(l) = j;
                    D.PSD.source.VIZU.F(l) = oV.F(k);
                    D.PSD.source.VIZU.labels{l} = table(k,1);
                    D.PSD.source.VIZU.callbacks(l) = oV.callbacks(k);
                end
            end
        end
    end
    if l >= 1
        D.other.val = l;
        D.PSD.source.VIZU.current = 1;
        D.PSD.source.VIZU.pst = unique(pst);
        D.PSD.source.VIZU.timeCourses = 1;
    else
        try,D.other = rmfield(D.other,'val');end
        D.PSD.source.VIZU.current = 0;
    end
end
set(D.PSD.handles.hfig,'userdata',D)
spm_eeg_review_callbacks('visu','main','info',D.PSD.VIZU.info)
drawnow
spm('pointer','arrow');

