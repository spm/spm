function [varargout] = spm_eeg_review_callbacks(varargin)
% Callbacks of the M/EEG Review facility
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jean Daunizeau
% $Id: spm_eeg_review_callbacks.m 3646 2009-12-16 16:01:03Z jean $

spm('pointer','watch');
drawnow expose

try
    D = get(spm_figure('FindWin','Graphics'),'userdata');
    handles = D.PSD.handles;
end



switch varargin{1}


    %% File I/O
    case 'file'
        switch varargin{2}
            case 'save'
                D0 = D;
                D = meeg(rmfield(D,'PSD'));
                save(D);
                D = D0;
                D.PSD.D0 = rmfield(D,'PSD');
                set(D.PSD.handles.hfig,'userdata',D)
                set(D.PSD.handles.BUTTONS.pop1,...
                    'BackgroundColor',[0.8314 0.8157 0.7843])
                
            case 'saveHistory'
                spm_eeg_history(D);
        end

    %% Get information from MEEG object
    case 'get'

        switch varargin{2}
            case 'VIZU'
                visuSensors             = varargin{3};
                VIZU.visuSensors        = visuSensors;
                VIZU.montage.clab       = {D.channels(visuSensors).label};
                if strcmp(D.transform.ID,'time')
                    M                       = sparse(length(visuSensors),length(D.channels));
                    M(sub2ind(size(M),1:length(visuSensors),visuSensors(:)')) = 1;
                    nts                     = min([2e2,D.Nsamples]);
                    decim                   = max([floor(D.Nsamples./nts),1]);
                    data                    = D.data.y(visuSensors,1:decim:D.Nsamples,:);
                    sd                      = mean(abs(data(:)-mean(data(:))));%std(data(:));
                    offset                  = (0:1:length(visuSensors)-1)'*(sd+eps)/2;
                    v_data                  = 0.25.*data +repmat(offset,[1 size(data,2) size(data,3)]);
                    ma                      = max(v_data(:))+sd;
                    mi                      = min(v_data(:))-sd;
                    ylim                    = [mi ma];
                    VIZU.visu_scale         = 0.25;
                    VIZU.FontSize           = 10;
                    VIZU.visu_offset        = sd;
                    VIZU.offset             = offset;
                    VIZU.ylim               = ylim;
                    VIZU.ylim0              = ylim;
                    VIZU.figname            = 'main visualization window';
                    VIZU.montage.M          = M;
                    VIZU.y2                 = permute(sum(data.^2,1),[2 3 1]);
                    VIZU.sci                = size(VIZU.y2,1)./D.Nsamples;
                else
                    nts                     = min([2e2,D.Nsamples*length(D.transform.frequencies)]);
                    decim                   = max([floor(D.Nsamples*length(D.transform.frequencies)./nts),1]);
                    data                    = D.data.y(visuSensors,:,1:decim:D.Nsamples,:);
                    VIZU.ylim               = [min(data(:)) max(data(:))];
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
            case 'uitable'
                D = getUItable(D);
            case 'prep'
                Finter = spm_figure('GetWin','Interactive');
                D = struct(get(Finter, 'UserData'));
                D0 = D.other(1).D0;
                D.other = rmfield(D.other,{'PSD','D0'});
                d1 = rmfield(D,'history');
                d0 = rmfield(D0,'history');
                if isequal(d1,d0)
                    % The objects only differ by their history
                    % => remove last operation from modified object
                    D.history(end) = [];
                end
                spm_eeg_review(D);
                hf = spm_figure('FindWin','Graphics');
                D = get(hf,'userdata');
                D.PSD.D0 = D0;
                set(hf,'userdata',D);
                spm_eeg_review_callbacks('visu','update')
                spm_clf(Finter)
        end


    %% Visualization callbacks
    case 'visu'

        switch varargin{2}

            %% Switch main uitabs: EEG/MEG/OTHER/INFO/SOURCE
            case 'main'

                try
                    D.PSD.VIZU.fromTab = D.PSD.VIZU.modality;
                catch
                    D.PSD.VIZU.fromTab = [];
                end

                switch varargin{3}
                    case 'eeg'
                        D.PSD.VIZU.modality = 'eeg';
                    case 'meg'
                        D.PSD.VIZU.modality = 'meg';
                    case 'megplanar'
                        D.PSD.VIZU.modality = 'megplanar';
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
                try
                    updateDisp(D,1);
                catch
                    set(D.PSD.handles.hfig,'userdata',D);
                end

            %% Switch from 'standard' to 'scalp' display type
            case 'switch'
                
                mod = get(gcbo,'userdata');
                if ~isequal(mod,D.PSD.VIZU.type)
                    if mod == 1
                        spm_eeg_review_callbacks('visu','main','standard')
                    else
                        spm_eeg_review_callbacks('visu','main','scalp')
                    end
                end

            %% Update display
            case 'update'

                try D = varargin{3};end
                updateDisp(D)

            %% Scalp interpolation
            case 'scalp_interp'

                if ~isempty([D.channels(:).X_plot2D])
                    x = round(mean(get(handles.axes(1),'xlim')));
                    ylim = get(handles.axes(1),'ylim');
                    if D.PSD.VIZU.type==1
                        in.hl = line('parent',handles.axes,...
                            'xdata',[x;x],...
                            'ydata',[ylim(1);ylim(2)]);
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
                        case 'megplanar'
                            I = D.PSD.MEGPLANAR.I;
                            in.type = 'MEGPLANAR';
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
                        y = y(:,x);
                        spm_eeg_plotScalpData(y,pos,labels,in);
                        try
                            D.PSD.handles.hli = in.hl;
                            set(D.PSD.handles.hfig,'userdata',D);
                        end
                    catch
                        msgbox('Get 2d positions for these channels!')
                    end
                else
                    msgbox('Get 2d positions for EEG/MEG channels!')
                end

            %% Display sensor positions (and canonical cortical mesh)
            case 'sensorPos'

                % get canonical mesh
                mco = fullfile(spm('Dir'),'canonical','cortex_5124.surf.gii');
                msc = fullfile(spm('Dir'),'canonical','scalp_2562.surf.gii');               
                
                % get and plot 3D sensor positions
                
                try     % EEG
                    try
                        for i=1:numel(D.other.inv{end}.datareg)
                            if isequal(D.other.inv{end}.datareg(i).modality,'EEG')
                                pos3d = spm_eeg_inv_transform_points(...
                                    D.other.inv{end}.datareg(i).toMNI,...
                                    D.other.inv{end}.datareg(i).sensors.pnt);
                            end
                        end
                        opt.figname = 'Coregistred EEG sensor positions';
                    catch
                        pos3d = [D.sensors.eeg.pnt];
                        pos3d = pos3d(D.PSD.EEG.I,:);
                        opt.figname = 'Uncoregistred EEG sensor positions';
                    end
                    pos3d(1,:);
                    % display canonical mesh
                    o = spm_eeg_render(mco,opt);
                    opt.hfig = o.handles.fi;
                    opt.ParentAxes = o.handles.ParentAxes;
                    o = spm_eeg_render(msc,opt);
                    set(o.handles.p,'FaceAlpha',0.75)
                    set(o.handles.transp,'value',0.75)
                    % display sensor position
                    figure(o.handles.fi);
                    set(opt.ParentAxes,'nextplot','add')
                    plot3(opt.ParentAxes,...
                        pos3d(:,1),pos3d(:,2),pos3d(:,3),'.');
                    try
                        labels = D.PSD.EEG.VIZU.montage.clab;
                        text(pos3d(:,1),pos3d(:,2),pos3d(:,3),...
                            labels,...
                            'parent',opt.ParentAxes);
                    end
                    axis(opt.ParentAxes,'equal')
                    axis(opt.ParentAxes,'tight')
                    axis(opt.ParentAxes,'off')
                end
                try     % MEG
                    clear opt pos3d o labels
                    try % multimodal EEG/MEG
                        for i=1:numel(D.other.inv{end}.datareg)
                            if isequal(D.other.inv{end}.datareg(i).modality,'MEG')
                                pos3d = spm_eeg_inv_transform_points(...
                                    D.other.inv{end}.datareg(i).toMNI,...
                                    D.other.inv{end}.datareg(i).sensors.pnt);
                            end
                        end
                        opt.figname = 'Coregistred MEG sensor positions';
                    catch
                        pos3d = [D.sensors.meg.pnt];
                        opt.figname = 'Uncoregistred MEG sensor positions';
                    end
                    pos3d(1,:);
                    % display canonical mesh
                    o = spm_eeg_render(mco,opt);
                    opt.hfig = o.handles.fi;
                    opt.ParentAxes = o.handles.ParentAxes;
                    o = spm_eeg_render(msc,opt);
                    set(o.handles.p,'FaceAlpha',0.75)
                    set(o.handles.transp,'value',0.75)
                    % display sensor position
                    figure(o.handles.fi);
                    set(opt.ParentAxes,'nextplot','add')
                    plot3(opt.ParentAxes,...
                        pos3d(:,1),pos3d(:,2),pos3d(:,3),'.');
                    try
                        labels = cat(2,...
                            D.PSD.MEG.VIZU.montage.clab,...
                            D.PSD.MEGPLANAR.VIZU.montage.clab);
                        text(pos3d(:,1),pos3d(:,2),pos3d(:,3),...
                            labels,...
                            'parent',opt.ParentAxes);
                    end
                    axis(opt.ParentAxes,'equal')
                    axis(opt.ParentAxes,'tight')
                    axis(opt.ParentAxes,'off')
                end
            
            %% Update display for 'SOURCE' main tab
            case 'inv'

                cla(D.PSD.handles.axes2,'reset')
                D.PSD.source.VIZU.current = varargin{3};
                updateDisp(D);

            %% Check xlim when resizing display window using 'standard'
            %% display type
            case 'checkXlim'

                xlim = varargin{3};
                ud = get(D.PSD.handles.gpa,'userdata');
                xm = mean(xlim);
                sw = abs(diff(xlim));
                if sw <= ud.v.minSizeWindow
                    sw = ud.v.minSizeWindow;
                elseif sw >= ud.v.nt
                    sw = ud.v.maxSizeWindow;
                elseif sw >= ud.v.maxSizeWindow
                    sw = ud.v.maxSizeWindow;
                end
                if xlim(1) <= 1 && xlim(end) >= ud.v.nt
                    xlim = [1,ud.v.nt];
                elseif xlim(1) <= 1
                    xlim = [1,sw];
                elseif xlim(end) >= ud.v.nt
                    xlim = [ud.v.nt-sw+1,ud.v.nt];
                end

                % Restrain buttons usage:
                if isequal(xlim,[1,ud.v.nt])
                    set(D.PSD.handles.BUTTONS.vb3,'enable','off')
                    set(handles.BUTTONS.slider_step,'visible','off')
                    set(D.PSD.handles.BUTTONS.goPlusOne,'visible','off');
                    set(D.PSD.handles.BUTTONS.goMinusOne,'visible','off');
                else
                    set(handles.BUTTONS.slider_step,...
                        'min',sw/2,'max',ud.v.nt-sw/2+1,...
                        'value',mean(xlim),...
                        'sliderstep',.1*[sw/(ud.v.nt-1) 4*sw/(ud.v.nt-1)],...
                        'visible','on');
                    set(D.PSD.handles.BUTTONS.goPlusOne,'visible','on');
                    set(D.PSD.handles.BUTTONS.goMinusOne,'visible','on');
                    if isequal(sw,ud.v.maxSizeWindow)
                        set(D.PSD.handles.BUTTONS.vb3,'enable','off')
                        set(D.PSD.handles.BUTTONS.vb4,'enable','on')
                    elseif isequal(sw,ud.v.minSizeWindow)
                        set(D.PSD.handles.BUTTONS.vb4,'enable','off')
                        set(D.PSD.handles.BUTTONS.vb3,'enable','on')
                    else
                        set(D.PSD.handles.BUTTONS.vb4,'enable','on')
                        set(D.PSD.handles.BUTTONS.vb3,'enable','on')
                    end
                    if xlim(1) == 1
                        set(D.PSD.handles.BUTTONS.goMinusOne,...
                            'visible','on','enable','off');
                        set(D.PSD.handles.BUTTONS.goPlusOne,...
                            'visible','on','enable','on');
                    elseif xlim(2) == ud.v.nt
                        set(D.PSD.handles.BUTTONS.goPlusOne,...
                            'visible','on','enable','off');
                        set(D.PSD.handles.BUTTONS.goMinusOne,...
                            'visible','on','enable','on');
                    else
                        set(D.PSD.handles.BUTTONS.goPlusOne,...
                            'visible','on','enable','on');
                        set(D.PSD.handles.BUTTONS.goMinusOne,...
                            'visible','on','enable','on');
                    end
                end
                if nargout >= 1
                    varargout{1} = xlim;
                else
                    D.PSD.VIZU.xlim = xlim;
                    set(D.PSD.handles.hfig,'userdata',D)
                end





            %% Contrast/intensity rescaling
            case 'iten_sc'

                switch D.PSD.VIZU.modality
                    case 'eeg'
                        D.PSD.EEG.VIZU.visu_scale = varargin{3}*D.PSD.EEG.VIZU.visu_scale;
                    case 'meg'
                        D.PSD.MEG.VIZU.visu_scale = varargin{3}*D.PSD.MEG.VIZU.visu_scale;
                    case 'megplanar'
                        D.PSD.MEGPLANAR.VIZU.visu_scale = varargin{3}*D.PSD.MEGPLANAR.VIZU.visu_scale;
                    case 'other'
                        D.PSD.other.VIZU.visu_scale = varargin{3}*D.PSD.other.VIZU.visu_scale;
                end
                updateDisp(D,3);


            %% Resize plotted data window ('standard' display type)
            case 'time_w'

                % Get current plotted data window range and limits
                xlim = get(handles.axes(1),'xlim');

                sw = varargin{3}*diff(xlim);
                xm = mean(xlim);
                xlim = xm + 0.5*[-sw,sw];

                xlim = spm_eeg_review_callbacks('visu','checkXlim',xlim);
                D.PSD.VIZU.xlim = xlim;

                updateDisp(D,4)


            %% Scroll through data ('standard' display type)
            case 'slider_t'

                offset = get(gco,'value');
                updateDisp(D)


            %% Scroll through data page by page  ('standard' display type)
            case 'goOne'

                % Get current plotted data window range and limits
                xlim = get(handles.axes(1),'xlim');
                sw = diff(xlim);
                xlim = xlim +varargin{3}*sw;
                xlim = spm_eeg_review_callbacks('visu','checkXlim',xlim);
                D.PSD.VIZU.xlim = xlim;
                updateDisp(D,4)


            %% Zoom
            case 'zoom'

                switch D.PSD.VIZU.type

                    case 1 % 'standard' display type

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

                    case 2 % 'scalp' display type

                        set(D.PSD.handles.BUTTONS.vb5,'value',1)
                        switch D.PSD.VIZU.modality
                            case 'eeg'
                                VIZU = D.PSD.EEG.VIZU;
                            case 'meg'
                                VIZU = D.PSD.MEG.VIZU;
                            case 'megplanar'
                                VIZU = D.PSD.MEGPLANAR.VIZU;
                            case 'other'
                                VIZU = D.PSD.other.VIZU;
                        end
                        try axes(D.PSD.handles.scale);end
                        [x] = ginput(1);
                        indAxes = get(gco,'userdata');
                        if ~~indAxes
                            hf = figure('color',[1 1 1]);
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
                                col = lines;
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
                                unit = 'unknown';
                                try
                                    unit = D.channels(VIZU.visuSensors(indAxes)).units;
                                end
                                if isequal(unit,'unknown')
                                    ylabel(ha2,'field intensity ')
                                else
                                    ylabel(ha2,['field intensity (in ',unit,')'])
                                end
                                title(ha2,['channel ',chanLabel,...
                                    ' (',D.channels(VIZU.visuSensors(indAxes)).type,')'])

                            else % time-frequency data

                                datai = squeeze(D.data.y(VIZU.visuSensors(indAxes),:,:,trN(1)));    
                                pst = (0:1/D.Fsample:(D.Nsamples-1)/D.Fsample) + D.timeOnset;
                                pst = pst*1e3;  % in msec
                                if any(size(datai)==1)
                                    hp2 = plot(datai,...
                                        'parent',ha2);
                                    set(ha2,'xtick',1:10:length(pst),'xticklabel',pst(1:10:length(pst)),...
                                        'xlim',[1 length(pst)]);
                                    xlabel(ha2,'time (in ms after time onset)')
                                    ylabel(ha2,'power in frequency space')
                                    title(ha2,['channel ',chanLabel,...
                                        ' (',D.channels(VIZU.visuSensors(indAxes)).type,')',...
                                        ' -- frequency: ',num2str(D.transform.frequencies),' Hz'])
                                else
                                    nx = max([1,length(pst)./10]);
                                    xtick = floor(1:nx:length(pst));
                                    ny = max([1,length(D.transform.frequencies)./10]);
                                    ytick = floor(1:ny:length(D.transform.frequencies));
                                    hp2 = image(datai,...
                                        'CDataMapping','scaled',...
                                        'parent',ha2);
                                    colormap(ha2,jet)
                                    colorbar('peer',ha2)
                                    set(ha2,...
                                        'xtick',xtick,...
                                        'xticklabel',pst(xtick),...
                                        'xlim',[0.5 length(pst)+0.5],...
                                        'ylim',[0.5 size(datai,1)+0.5],...
                                        'ytick',ytick,...
                                        'yticklabel',D.transform.frequencies(ytick));
                                    xlabel(ha2,'time (in ms after time onset)')
                                    ylabel(ha2,'frequency (in Hz)')
                                    title(ha2,['channel ',chanLabel,...
                                        ' (',D.channels(VIZU.visuSensors(indAxes)).type,')'])
                                    caxis(ha2,VIZU.ylim)
                                end

                            end

                            axes(ha2)
                        end
                        set(D.PSD.handles.BUTTONS.vb5,'value',0)
                end

            otherwise;disp('unknown command !')


        end

    %% Events callbacks accessible from uicontextmenu
    %% ('standard' display type when playing with 'continuous' data)
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

                    updateDisp(D,2,currentEvent)

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
                        updateDisp(D,4)
                    end
                end



            % Execute actions accessible from the event contextmenu : delete event
            case 'deleteEvent'

                D.trials.events(currentEvent) = [];
                updateDisp(D,2)

        end

    %% Events callbacks
    case 'select'

        switch varargin{2}

            %% Switch to another trial (when playing with 'epoched' data)
            case 'switch'
                trN = get(gco,'value');
                if ~strcmp(D.PSD.VIZU.modality,'source') && D.PSD.VIZU.type == 2
                    handles = rmfield(D.PSD.handles,'PLOT');
                    D.PSD.handles = handles;
                else
                    try cla(D.PSD.handles.axes2,'reset');end
                end
                D.PSD.trials.current = trN;
                status = any([D.trials(trN).bad]);
                try
                    if status
                        str = 'declare as not bad';
                    else
                        str = 'declare as bad';
                    end
                    ud = get(D.PSD.handles.BUTTONS.badEvent,'userdata');
                    set(D.PSD.handles.BUTTONS.badEvent,...
                        'tooltipstring',str,...
                        'cdata',ud.img{2-status},'userdata',ud)
                end
                updateDisp(D,1)

            %% Switch event to 'bad' (when playing with 'epoched' data)
            case 'bad'
                trN = D.PSD.trials.current;
                status = any([D.trials(trN).bad]);
                str1 = 'not bad';
                str2 = 'bad';
                if status
                    bad = 0;
                    lab = [' (',str1,')'];
                    str = ['declare as ',str2];
                else
                    bad = 1;
                    lab = [' (',str2,')'];
                    str = ['declare as ',str1];
                end
                nt = length(trN);
                for i=1:nt
                    D.trials(trN(i)).bad = bad;
                    D.PSD.trials.TrLabels{trN(i)} = ['Trial ',num2str(trN(i)),...
                        ': ',D.trials(trN(i)).label,lab];
                end
                set(D.PSD.handles.BUTTONS.list1,'string',D.PSD.trials.TrLabels);
                ud = get(D.PSD.handles.BUTTONS.badEvent,'userdata');
                set(D.PSD.handles.BUTTONS.badEvent,...
                    'tooltipstring',str,...
                    'cdata',ud.img{2-bad},'userdata',ud)
                set(D.PSD.handles.hfig,'userdata',D)

            %% Add an event to current selection
            %% (when playing with 'continuous' data)
            case 'add'
                [x,tmp] = ginput(1);
                x = round(x);
                x(1) = min([max([1 x(1)]) D.Nsamples]);
                Nevents = length(D.trials.events);
                D.trials.events(Nevents+1).time = x./D.Fsample;
                D.trials.events(Nevents+1).duration = 0;
                D.trials.events(Nevents+1).type = 'Manual';
                D.PSD.handles.PLOT.e(Nevents+1) = 0;
                if Nevents > 0
                    D.trials.events(Nevents+1).value = D.trials.events(Nevents).value;
                else
                    D.trials.events(Nevents+1).value = 0;
                end
                % Enable tools on selections
                set(handles.BUTTONS.sb2,'enable','on');
                set(handles.BUTTONS.sb3,'enable','on');
                % Update display
                updateDisp(D,2,Nevents+1)


            %% scroll through data upto next event
            %% (when playing with 'continuous' data)
            case 'goto'
                here                    = get(handles.BUTTONS.slider_step,'value');
                x                       = [D.trials.events.time]';
                xm                      = x.*D.Fsample;
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
                        updateDisp(D,4)
                    end
                end

        end

   %% Edit callbacks (from spm_eeg_prep_ui)
    case 'edit'

        switch varargin{2}

            case 'prep'

                try rotate3d off;end
                spm_eeg_prep_ui;
                Finter = spm_figure('GetWin','Interactive');
                D0 = D.PSD.D0;
                D = rmfield(D,'PSD');
                if isempty(D.other)
                    D.other = struct([]);
                end
                D.other(1).PSD = 1;
                D.other(1).D0 = D0;
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


% Check changes in the meeg object
if isstruct(D)&& isfield(D,'PSD') && ...
        isfield(D.PSD,'D0')
    d1 = rmfield(D,{'history','PSD'});
    d0 = rmfield(D.PSD.D0,'history');
    if isequal(d1,d0)
        set(D.PSD.handles.BUTTONS.pop1,...
            'BackgroundColor',[0.8314 0.8157 0.7843])
    else
        set(D.PSD.handles.BUTTONS.pop1,...
            'BackgroundColor',[1 0.5 0.5])
    end
end
spm('pointer','arrow');
drawnow expose

%% Main update display
function [] = updateDisp(D,flags,in)
% This function updates the display of the data and events.

if ~exist('flag','var')
    flag = 0;
end
if ~exist('in','var')
    in = [];
end
handles = D.PSD.handles;



% Create intermediary display variables : events
figure(handles.hfig)


% Get current event
try
    trN = D.PSD.trials.current;
catch
    trN = 1;
end

if ~strcmp(D.PSD.VIZU.modality,'source')

    switch D.PSD.VIZU.modality
        case 'eeg'
            VIZU = D.PSD.EEG.VIZU;
        case 'meg'
            VIZU = D.PSD.MEG.VIZU;
        case 'megplanar'
            VIZU = D.PSD.MEGPLANAR.VIZU;
        case 'other'
            VIZU = D.PSD.other.VIZU;
        case 'info'
            return
    end


    switch D.PSD.VIZU.type

        case 1

            % Create new data to display
            %   - switch from scalp to standard displays
            %   - switch from EEG/MEG/OTHER/info/inv
            if ismember(1,flags)
                % delete previous axes...
                try
                    delete(D.PSD.handles.axes)
                    delete(D.PSD.handles.gpa)
                    delete(D.PSD.handles.BUTTONS.slider_step)
                end
                % gather info for core display function
                options.hp = handles.hfig;
                options.Fsample = D.Fsample;
                options.timeOnset = D.timeOnset;
                options.M = VIZU.visu_scale*full(VIZU.montage.M);
                options.bad = [D.channels(VIZU.visuSensors(:)).bad];
                if strcmp(D.PSD.type,'continuous') && ~isempty(D.trials.events)
                    trN = 1;
                    Nevents = length(D.trials.events);
                    x1 = {D.trials.events(:).type}';
                    x2 = {D.trials.events(:).value}';
                    if ~iscellstr(x1)
                        [y1,i1,j1] = unique(cell2mat(x1));
                    else
                        [y1,i1,j1] = unique(x1);
                    end
                    if ~iscellstr(x2)
                        [y2,i2,j2] = unique(cell2mat(x2));
                    else
                        [y2,i2,j2] = unique(x2);
                    end
                    A = [j1(:),j2(:)];
                    [ya,ia,ja] = unique(A,'rows');
                    options.events = rmfield(D.trials.events,{'duration','value'});
                    for i=1:length(options.events)
                        options.events(i).time = options.events(i).time.*D.Fsample;% +1;
                        options.events(i).type = ja(i);
                    end
                end
                if strcmp(D.PSD.type,'continuous')
                    options.minSizeWindow = 200;
                    try
                        options.itw = round(D.PSD.VIZU.xlim(1):D.PSD.VIZU.xlim(2));
                    end
                elseif strcmp(D.PSD.type,'epoched')
                    options.minSizeWindow = 20;
                    try
                        options.itw = round(D.PSD.VIZU.xlim(1):D.PSD.VIZU.xlim(2));
                    catch
                        options.itw = 1:D.Nsamples;
                    end
                else
                    try
                        options.itw = round(D.PSD.VIZU.xlim(1):D.PSD.VIZU.xlim(2));
                    catch
                        options.itw = 1:D.Nsamples;
                    end
                    options.minSizeWindow = 20;
                end
                options.minY = min(VIZU.ylim)-eps;
                options.maxY = max(VIZU.ylim)+eps;
                options.ds = 5e2;
                options.pos1 = [0.08 0.11 0.86 0.79];
                options.pos2 = [0.08 0.07 0.86 0.025];
                options.pos3 = [0.08 0.02 0.86 0.02];
                options.maxSizeWindow = 1e5;
                options.tag = 'plotEEG';
                options.offset = VIZU.offset;
                options.ytick = VIZU.offset;
                options.yticklabel = VIZU.montage.clab;
                options.callback = ['spm_eeg_review_callbacks(''visu'',''checkXlim''',...
                    ',get(ud.v.handles.axes,''xlim''))'];
                % Use file_array for 'continuous' data.
                if strcmp(D.PSD.type,'continuous')
                    options.transpose = 1;
                    ud = spm_DisplayTimeSeries(D.data.y,options);
                else
                    ud = spm_DisplayTimeSeries(D.data.y(:,:,trN(1))',options);
                end
                % update D
                D.PSD.handles.axes = ud.v.handles.axes;
                D.PSD.handles.gpa = ud.v.handles.gpa;
                D.PSD.handles.BUTTONS.slider_step = ud.v.handles.hslider;
                D.PSD.handles.PLOT.p = ud.v.handles.hp;
                % Create uicontextmenu for events (if any)
                if isfield(options,'events')
                    D.PSD.handles.PLOT.e = [ud.v.et(:).hp];
                    axes(D.PSD.handles.axes)
                    for i=1:length(options.events)
                        sc.currentEvent = i;
                        sc.eventType    = D.trials(trN(1)).events(i).type;
                        sc.eventValue   = D.trials(trN(1)).events(i).value;
                        sc.N_select     = Nevents;
                        psd_defineMenuEvent(D.PSD.handles.PLOT.e(i),sc);
                    end
                end
                for i=1:length(D.PSD.handles.PLOT.p)
                    cmenu = uicontextmenu;
                    uimenu(cmenu,'Label',['channel ',num2str(VIZU.visuSensors(i)),': ',VIZU.montage.clab{i}]);
                    uimenu(cmenu,'Label',['type: ',D.channels(VIZU.visuSensors(i)).type]);
                    uimenu(cmenu,'Label',['bad: ',num2str(D.channels(VIZU.visuSensors(i)).bad)],...
                        'callback',@switchBC,'userdata',i,...
                        'BusyAction','cancel',...
                        'Interruptible','off');
                    set(D.PSD.handles.PLOT.p(i),'uicontextmenu',cmenu);
                end
                set(D.PSD.handles.hfig,'userdata',D);
                spm_eeg_review_callbacks('visu','checkXlim',...
                    get(D.PSD.handles.axes,'xlim'))
            end

            % modify events properties (delete,add,time,...)
            if ismember(2,flags)
                Nevents = length(D.trials.events);
                if Nevents < length(D.PSD.handles.PLOT.e)
                    action = 'delete';
                    try,delete(D.PSD.handles.PLOT.e),end
                    try,D.PSD.handles.PLOT.e = [];end
                else
                    action = 'modify';
                end
                col = lines;
                col = col(1:7,:);
                x1 = {D.trials.events(:).type}';
                x2 = {D.trials.events(:).value}';
                if ~iscellstr(x1)
                    [y1,i1,j1] = unique(cell2mat(x1));
                else
                    [y1,i1,j1] = unique(x1);
                end
                if ~iscellstr(x2)
                    [y2,i2,j2] = unique(cell2mat(x2));
                else
                    [y2,i2,j2] = unique(x2);
                end
                A = [j1(:),j2(:)];
                [ya,ia,ja] = unique(A,'rows');
                events = rmfield(D.trials.events,{'duration','value'});
                switch action
                    case 'delete'
                        %spm_progress_bar('Init',Nevents,'Replacing events');
                        axes(D.PSD.handles.axes)
                        for i=1:Nevents
                            events(i).time = D.trials.events(i).time.*D.Fsample;% +1;
                            events(i).type = ja(i);
                            events(i).col = mod(events(i).type+7,7)+1;
                            D.PSD.handles.PLOT.e(i) = plot(D.PSD.handles.axes,...
                                events(i).time.*[1 1],...
                                VIZU.ylim,...
                                'color',col(events(i).col,:),...
                                'userdata',i,...
                                'ButtonDownFcn','set(gco,''selected'',''on'')',...
                                'Clipping','on');
                            % Add events uicontextmenu
                            sc.currentEvent = i;
                            sc.eventType    = D.trials(trN(1)).events(i).type;
                            sc.eventValue   = D.trials(trN(1)).events(i).value;
                            sc.N_select     = Nevents;
                            psd_defineMenuEvent(D.PSD.handles.PLOT.e(i),sc);
                            %spm_progress_bar('Set',i)
                        end
                        %spm_progress_bar('Clear')
                    case 'modify'
                        events(in).time = D.trials.events(in).time.*D.Fsample;% +1;
                        events(in).type = ja(in);
                        events(in).col = mod(events(in).type+7,7)+1;
                        D.PSD.handles.PLOT.e(in) = plot(D.PSD.handles.axes,events(in).time.*[1 1],...
                            VIZU.ylim,'color',col(events(in).col,:));
                        set(D.PSD.handles.PLOT.e(in),'userdata',in,...
                            'ButtonDownFcn','set(gco,''selected'',''on'')',...
                            'Clipping','on');
                        % Add events uicontextmenu
                        sc.currentEvent = in;
                        sc.eventType    = D.trials(trN(1)).events(in).type;
                        sc.eventValue   = D.trials(trN(1)).events(in).value;
                        sc.N_select     = Nevents;
                        psd_defineMenuEvent(D.PSD.handles.PLOT.e(in),sc);
                end     
                set(handles.hfig,'userdata',D);
            end

            % modify scaling factor
            if ismember(3,flags)
                ud = get(D.PSD.handles.gpa,'userdata');
                ud.v.M = VIZU.visu_scale*full(VIZU.montage.M);
                xw = floor(get(ud.v.handles.axes,'xlim'));
                xw(1) = max([1,xw(1)]);
                if ~ud.v.transpose
                    My = ud.v.M*ud.y(xw(1):1:xw(2),:)';
                else
                    My = ud.v.M*ud.y(:,xw(1):1:xw(2));
                end
                for i=1:ud.v.nc
                    set(ud.v.handles.hp(i),'xdata',xw(1):1:xw(2),'ydata',My(i,:)+ud.v.offset(i))
                end
                set(ud.v.handles.axes,'ylim',[ud.v.mi ud.v.ma],'xlim',xw);
                set(D.PSD.handles.gpa,'userdata',ud);
                set(handles.hfig,'userdata',D);
            end

            % modify plotted time window (goto, ...)
            if ismember(4,flags)
                ud = get(D.PSD.handles.gpa,'userdata');
                xw = floor(D.PSD.VIZU.xlim);
                xw(1) = max([1,xw(1)]);
                if ~ud.v.transpose
                    My = ud.v.M*ud.y(xw(1):1:xw(2),:)';
                else
                    My = ud.v.M*ud.y(:,xw(1):1:xw(2));
                end
                for i=1:ud.v.nc
                    set(ud.v.handles.hp(i),'xdata',xw(1):1:xw(2),'ydata',My(i,:)+ud.v.offset(i))
                end
                set(ud.v.handles.axes,'ylim',[ud.v.mi ud.v.ma],'xlim',xw);
                set(ud.v.handles.pa,'xdata',[xw,fliplr(xw)]);
                set(ud.v.handles.lb,'xdata',[xw(1) xw(1)]);
                set(ud.v.handles.rb,'xdata',[xw(2) xw(2)]);
                sw = diff(xw);
                set(ud.v.handles.hslider,'value',mean(xw),...
                    'min',1+sw/2,'max',ud.v.nt-sw/2,...
                    'sliderstep',.1*[sw/(ud.v.nt-1) 4*sw/(ud.v.nt-1)]);
                set(handles.hfig,'userdata',D);
            end


        case 2

            if strcmp(D.transform.ID,'time')

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
                    
                    if miY == 0 && maY == 0
                        miY = -eps;
                        maY = eps;
                    else
                        miY = miY - miY.*1e-3;
                        maY = maY + maY.*1e-3;
                    end
                    
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
                    datai = squeeze(D.data.y(VIZU.visuSensors(i),:,:,trN(1)));
                    miY = min(datai(:));
                    maY = max(datai(:));
                    if any(size(datai)==1)
                        D.PSD.handles.PLOT.im(i) = plot(datai,...
                            'parent',handles.axes(i),...
                            'tag','plotEEG',...
                            'userdata',i,...
                            'hittest','off');
                        set(handles.axes(i),...
                            'ylim',[miY maY]);
                    else
                        D.PSD.handles.PLOT.im(i) = image(datai,...
                            'parent',handles.axes(i),...
                            'CDataMapping','scaled',...
                            'tag','plotEEG',...
                            'userdata',i,...
                            'hittest','off');
                    end
                    set(handles.fra(i),'uicontextmenu',cmenu);
                end
                colormap(jet)
                % This normalizes colorbars across channels and trials:
                for i=1:length(VIZU.visuSensors)
                    caxis(handles.axes(i),VIZU.ylim);
                end
                set(handles.hfig,'userdata',D);

            end
    end


else  % source space

    % get model/trial info
    VIZU = D.PSD.source.VIZU;
    isInv = VIZU.isInv;
    Ninv = length(isInv);
    invN = VIZU.isInv(D.PSD.source.VIZU.current);
    F  = VIZU.F;
    ID = VIZU.ID;
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
                  
        if Ninv>1
            if isnan(ID(invN))
                xF = find(isnan(ID));
            else
                xF = find(abs(ID-ID(invN))<eps);
            end
            if length(xF)>1
                D.PSD.handles.hbar = bar(D.PSD.handles.BMCplot,...
                    xF ,F(xF)-min(F(xF)),...
                    'barwidth',0.5,...
                    'FaceColor',0.5*[1 1 1],...
                    'visible','off',...
                    'tag','plotEEG');
                D.PSD.handles.BMCcurrent = plot(D.PSD.handles.BMCplot,...
                    find(xF==invN),0,'ro',...
                    'visible','off',...
                    'tag','plotEEG');
                set(D.PSD.handles.BMCplot,...
                    'xtick',xF,...
                    'xticklabel',D.PSD.source.VIZU.labels(xF),...
                    'xlim',[0,length(xF)+1]);
                drawnow
            else
                cla(D.PSD.handles.BMCplot);
                set(D.PSD.handles.BMCplot,...
                    'xtick',[],...
                    'xticklabel',{});
            end
        end
        
        % get model/trial time series
        D.PSD.source.VIZU.J = zeros(model.Nd,size(model.T,1));
        D.PSD.source.VIZU.J(model.Is,:) = model.J{trN(1)}*model.T';
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
                Jp(1,:) = min(D.PSD.source.VIZU.J,[],1);
                Jp(2,:) = max(D.PSD.source.VIZU.J,[],1);
                D.PSD.source.VIZU.plotTC = plot(D.PSD.handles.axes2,...
                    model.pst,Jp','color',0.5*[1 1 1]);
                set(D.PSD.handles.axes2,'hittest','off')
                % Add virtual electrode
                %                 try
                %                     ve = D.PSD.source.VIZU.ve;
                %                 catch
                [mj ve] = max(max(abs(D.PSD.source.VIZU.J),[],2));
                D.PSD.source.VIZU.ve =ve;
                %                 end
                Jve = D.PSD.source.VIZU.J(D.PSD.source.VIZU.ve,:);
                set(D.PSD.handles.axes2,'nextplot','add')
                try
                    qC  = model.qC(ve).*diag(model.qV)';
                    ci  = 1.64*sqrt(qC);
                    D.PSD.source.VIZU.pve2 = plot(D.PSD.handles.axes2,...
                        model.pst,Jve +ci,'b:',model.pst,Jve -ci,'b:');
                end
                D.PSD.source.VIZU.pve = plot(D.PSD.handles.axes2,...
                    model.pst,Jve,'color','b');
                set(D.PSD.handles.axes2,'nextplot','replace')
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
    case 'megplanar'
        I = D.PSD.MEGPLANAR.I;
        VIZU = D.PSD.MEGPLANAR.VIZU;    
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
        %         ud = get(D.PSD.handles.axes);
        %         ud.v.bad(ind) = status;
        %         set(D.PSD.handles.axes,'userdata',ud);
    case 2
        set(D.PSD.handles.axes(ind),'Color',color);
        set(D.PSD.handles.fra(ind),'uicontextmenu',cmenu);
        set(D.PSD.handles.PLOT.p(:,ind),'uicontextmenu',cmenu);
        axes(D.PSD.handles.scale)
end

d1 = rmfield(D,{'history','PSD'});
d0 = rmfield(D.PSD.D0,'history');
if isequal(d1,d0)
    set(D.PSD.handles.BUTTONS.pop1,...
        'BackgroundColor',[0.8314 0.8157 0.7843])
else
    set(D.PSD.handles.BUTTONS.pop1,...
        'BackgroundColor',[1 0.5 0.5])
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
if isfield(D.other.inv{invN}.inverse, 'modality')
    mod0 = D.other.inv{invN}.inverse.modality;
    if ischar(mod0)
        mod = mod0;
    else
        mod = [];
        for i = 1:length(mod0)
            mod = [mod,' ',mod0{i}];
        end
    end
    str{3} = ['Modality: ',mod];
else % For backward compatibility
    try
        mod0 = D.other.inv{invN}.modality;
        if ischar(mod0)
            mod = mod0;
        else
            mod = [];
            for i = 1:length(mod0)
                mod = [mod,' ',mod0{i}];
            end
        end
        str{3} = ['Modality: ',mod];
    catch
        str{3} = 'Modality: ?';
    end
end
    
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
str{1} = ['File name: ',fullfile(D.path,D.fname)];
str{2} = ['Type: ',D.type];
if ~strcmp(D.transform.ID,'time')
    str{2} = [str{2},' (time-frequency data, from ',...
        num2str(D.transform.frequencies(1)),'Hz to ',...
        num2str(D.transform.frequencies(end)),'Hz'];
    if strcmp(D.transform.ID,'TF')
        str{2} = [str{2},')'];
    else
        str{2} = [str{2},': phase)'];
    end
end
delta_t = D.Nsamples./D.Fsample;
gridTime = (1:D.Nsamples)./D.Fsample + D.timeOnset;
str{3} = ['Number of time samples: ',num2str(D.Nsamples),' (',num2str(delta_t),' sec, from ',...
    num2str(gridTime(1)),'s to ',num2str(gridTime(end)),'s)'];
str{4} = ['Time sampling frequency: ',num2str(D.Fsample),' Hz'];
nb = length(find([D.channels.bad]));
str{5} = ['Number of channels: ',num2str(length(D.channels)),' (',num2str(nb),' bad channels)'];
nb = length(find([D.trials.bad]));
if strcmp(D.type,'continuous')
    if isfield(D.trials(1),'events')
        str{6} = ['Number of events: ',num2str(length(D.trials(1).events))];
    else
        str{6} = ['Number of events: ',num2str(0)];
    end
else
    str{6} = ['Number of trials: ',num2str(length(D.trials)),' (',num2str(nb),' bad trials)'];
end
% try,str{7} = ['Time onset: ',num2str(D.timeOnset),' sec'];end




%% extracting data from spm_uitable java object
function [D] = getUItable(D)
ht = D.PSD.handles.infoUItable;
cn = get(ht,'columnNames');
table = get(ht,'data');
% !! there is some redundancy here --> to be optimized...
table2 = spm_uitable('get',ht);
emptyTable = 0;
try
    emptyTable = isempty(cell2mat(table2));
end


if length(cn) == 5  % channel info
    if ~emptyTable
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
                    case 'megplanar'
                        D.channels(i).type = 'MEGPLANAR';
                    case 'megmag'
                        D.channels(i).type = 'MEGMAG';
                    case 'meggrad'
                        D.channels(i).type = 'MEGGRAD';
                    case 'refmag'
                        D.channels(i).type = 'REFMAG';
                    case 'refgrad'
                        D.channels(i).type = 'REFGRAD';     
                    case 'lfp'
                        D.channels(i).type = 'LFP';
                    case 'eog'
                        D.channels(i).type = 'EOG';
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
        D.PSD.MEG.I  = sort([find(strcmp('MEGMAG',{D.channels.type})),...
            find(strcmp('MEGGRAD',{D.channels.type})) find(strcmp('MEG',{D.channels.type}))]);
        D.PSD.MEGPLANAR.I  = find(strcmp('MEGPLANAR',{D.channels.type}));
        D.PSD.other.I = setdiff(1:nc,[D.PSD.EEG.I(:);D.PSD.MEG.I(:)]);
        if ~isempty(D.PSD.EEG.I)
            [out] = spm_eeg_review_callbacks('get','VIZU',D.PSD.EEG.I);
            D.PSD.EEG.VIZU = out;
        else
            D.PSD.EEG.VIZU = [];
        end
        if ~isempty(D.PSD.MEG.I)
            [out] = spm_eeg_review_callbacks('get','VIZU',D.PSD.MEG.I);
            D.PSD.MEG.VIZU = out;
        else
            D.PSD.MEG.VIZU = [];
        end
        if ~isempty(D.PSD.MEGPLANAR.I)
            [out] = spm_eeg_review_callbacks('get','VIZU',D.PSD.MEGPLANAR.I);
            D.PSD.MEGPLANAR.VIZU = out;
        else
            D.PSD.MEGPLANAR.VIZU = [];
        end
        if ~isempty(D.PSD.other.I)
            [out] = spm_eeg_review_callbacks('get','VIZU',D.PSD.other.I);
            D.PSD.other.VIZU = out;
        else
            D.PSD.other.VIZU = [];
        end
    else

    end
elseif length(cn) == 7
    if strcmp(D.type,'continuous')
        if ~emptyTable
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
            D.trials(1).events = [];
            delete(ht);
        end
    else
        if ~emptyTable
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
                        D.trials(i).events.value = table(i,3);%str2double(table(i,3));
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
        else
        end
    end

elseif length(cn) == 3
    if ~emptyTable
        nt = length(D.trials);
        for i=1:nt
            if ~isempty(table(i,1))
                D.trials(i).label = table(i,1);
            end
            D.PSD.trials.TrLabels{i} = ['Trial ',num2str(i),' (average of ',...
                num2str(D.trials(i).repl),' events): ',D.trials(i).label];
        end
    else
    end

elseif length(cn) == 12     % source reconstructions

    if ~emptyTable
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
            try D.other = rmfield(D.other,'val');end
            D.PSD.source.VIZU.current = 0;
        end
    else
        try D.other = rmfield(D.other,'val');end
        try D.other = rmfield(D.other,'inv');end
        D.PSD.source.VIZU.current = 0;
        D.PSD.source.VIZU.isInv = [];
        D.PSD.source.VIZU.pst = [];
        D.PSD.source.VIZU.F = [];
        D.PSD.source.VIZU.labels = [];
        D.PSD.source.VIZU.callbacks = [];
        D.PSD.source.VIZU.timeCourses = [];
        delete(ht)
    end
end
set(D.PSD.handles.hfig,'userdata',D)
spm_eeg_review_callbacks('visu','main','info',D.PSD.VIZU.info)
