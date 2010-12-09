function [D] = spm_eeg_review_switchDisplay(D)
% Switch between displays in the M/EEG Review facility
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jean Daunizeau
% $Id: spm_eeg_review_switchDisplay.m 4136 2010-12-09 22:22:28Z guillaume $

try % only if already displayed stuffs
    handles = rmfield(D.PSD.handles,'PLOT');
    D.PSD.handles = handles;
end


switch D.PSD.VIZU.modality

    case 'source'
        delete(findobj('tag','plotEEG'));
        [D] = visuRecon(D);

    case 'info'
        
        [D] = DataInfo(D);
        set(D.PSD.handles.hfig,'userdata',D)

    otherwise % plot data (EEG/MEG/OTHER)

        try
            y = D.data.y(:,D.PSD.VIZU.xlim(1):D.PSD.VIZU.xlim(2));
            % ! accelerates memory mapping reading
        catch
            D.PSD.VIZU.xlim = [1,min([5e2,D.Nsamples])];
        end
        switch  D.PSD.VIZU.type

            case 1
                delete(findobj('tag','plotEEG'))
                [D] = standardData(D);
                cameratoolbar('resetcamera')
                try cameratoolbar('close'); end

            case 2
                delete(findobj('tag','plotEEG'))
                [D] = scalpData(D);
                cameratoolbar('resetcamera')
                try cameratoolbar('close'); end

        end

end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Standard EEG/MEG data plot
function [D] = standardData(D)

% POS = get(D.PSD.handles.hfig,'position');

switch D.PSD.VIZU.modality
    case 'eeg'
        I = D.PSD.EEG.I;
        scb = 6;
    case 'meg'
        I = D.PSD.MEG.I;
        scb = 6;
    case 'megplanar'
        I = D.PSD.MEGPLANAR.I;
        scb = 6;    
    case 'other'
        I = D.PSD.other.I;
        scb = [];  % no scalp interpolation button
end

if isempty(I)

    uicontrol('style','text',...
        'units','normalized','Position',[0.14 0.84 0.7 0.04],...
        'string','No channel of this type in the SPM data file !',...
        'BackgroundColor',0.95*[1 1 1],...
        'tag','plotEEG')

else

    if ~strcmp(D.transform.ID,'time')

        uicontrol('style','text',...
            'units','normalized','Position',[0.14 0.84 0.7 0.04],...
            'string','Not for time-frequency data !',...
            'BackgroundColor',0.95*[1 1 1],...
            'tag','plotEEG')

    else


        D.PSD.VIZU.type = 1;

        % add buttons
        object.type = 'buttons';
        object.options.multSelect = 0;
        object.list = [2;3;4;5;scb];
        switch D.PSD.type
            case 'continuous'
                object.list = [object.list;9];
            case 'epoched'
                object.list = [object.list;7;11];
                if strcmp(D.type,'single')
                    object.list = [object.list;13];
                end
        end
        D = spm_eeg_review_uis(D,object);


    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 'SPM-like' EEG/MEG data plot
function [D] = scalpData(D)

% POS = get(D.PSD.handles.hfig,'position');

switch D.PSD.VIZU.modality
    case 'eeg'
        I = D.PSD.EEG.I;
    case 'meg'
        I = D.PSD.MEG.I;
    case 'megplanar'
        I = D.PSD.MEGPLANAR.I;    
    case 'other'
        I = D.PSD.other.I;
end

if isempty(I)

    uicontrol('style','text',...
        'units','normalized','Position',[0.14 0.84 0.7 0.04],...
        'string','No channel of this type in the SPM data file !',...
        'BackgroundColor',0.95*[1 1 1],...
        'tag','plotEEG')

else

    if strcmp(D.PSD.type,'continuous')

        uicontrol('style','text',...
            'units','normalized','Position',[0.14 0.84 0.7 0.04],...
            'string','Only for epoched data !',...
            'BackgroundColor',0.95*[1 1 1],...
            'tag','plotEEG')

    else

        D.PSD.VIZU.type = 2;
        % add buttons
        object.type = 'buttons';
        object.list = [5;7];
        if strcmp(D.transform.ID,'time') % only for time data!
            object.options.multSelect = 1;
            object.list = [object.list;4;6;11];
        else
            object.options.multSelect = 0;
        end
        if strcmp(D.type,'single')
            object.list = [object.list;13];
        end
        D = spm_eeg_review_uis(D,object);

        % add axes (!!give channels!!)
        switch D.PSD.VIZU.modality
            case 'eeg'
                I = D.PSD.EEG.I;
                ylim = D.PSD.EEG.VIZU.ylim;
            case 'meg'
                I = D.PSD.MEG.I;
                ylim = D.PSD.MEG.VIZU.ylim;
            case 'megplanar'
                I = D.PSD.MEGPLANAR.I;
                ylim = D.PSD.MEGPLANAR.VIZU.ylim;
            case 'other'
                I = D.PSD.other.I;
                ylim = D.PSD.other.VIZU.ylim;
        end
        object.type = 'axes';
        object.what = 'scalp';
        object.options.channelPlot = I;
        object.options.ylim = ylim;
        D = spm_eeg_review_uis(D,object);

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% RENDERING OF INVERSE SOLUTIONS
function [D] = visuRecon(D)
POS = get(D.PSD.handles.hfig,'position');

if ~~D.PSD.source.VIZU.current
    
    isInv = D.PSD.source.VIZU.isInv;
    Ninv = length(isInv);
    if D.PSD.source.VIZU.current > Ninv
        D.PSD.source.VIZU.current = 1;
    end
    invN = isInv(D.PSD.source.VIZU.current);
    pst = D.PSD.source.VIZU.pst;
    F  = D.PSD.source.VIZU.F;
    ID = D.PSD.source.VIZU.ID;

    % create uitabs for inverse solutions
    hInv = D.PSD.handles.tabs.hp;
    [h] = spm_uitab(hInv,D.PSD.source.VIZU.labels,...
        D.PSD.source.VIZU.callbacks,'plotEEG',...
        D.PSD.source.VIZU.current);
    D.PSD.handles.SubTabs_inv = h;

    trN = D.PSD.trials.current(1);
    model = D.other.inv{invN}.inverse;
    D.PSD.source.VIZU.J = zeros(model.Nd,size(model.T,1));
    D.PSD.source.VIZU.J(model.Is,:) = model.J{trN}*model.T';
    D.PSD.source.VIZU.miJ = min(min(D.PSD.source.VIZU.J));
    D.PSD.source.VIZU.maJ = max(max(D.PSD.source.VIZU.J));

    J = D.PSD.source.VIZU.J;
    miJ = D.PSD.source.VIZU.miJ;
    maJ = D.PSD.source.VIZU.maJ;
    time = (model.pst-0).^2;
    indTime = find(time==min(time));
    gridTime = model.pst(indTime);

    % create axes
    object.type = 'axes';
    object.what = 'source';
    object.options.Ninv = Ninv;
    object.options.miJ = miJ;
    object.options.maJ = maJ;
    object.options.pst = pst;
    D = spm_eeg_review_uis(D,object);

    % plot BMC free energies in appropriate axes
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
        end
    end

    % Create mesh and related objects
    Dmesh = D.other.inv{invN}.mesh;
    mesh.vertices = Dmesh.tess_mni.vert;
    mesh.faces = Dmesh.tess_mni.face;
    options.texture = J(:,indTime);
    options.hfig = D.PSD.handles.hfig;
    options.ParentAxes = D.PSD.handles.axes;
    options.tag = 'plotEEG';
    options.visible = 'off';
    [out] = spm_eeg_render(mesh,options);
    D.PSD.handles.mesh = out.handles.p;
    D.PSD.handles.BUTTONS.transp = out.handles.transp;
    D.PSD.handles.colorbar = out.handles.hc;
    D.PSD.handles.BUTTONS.ct1 = out.handles.s1;
    D.PSD.handles.BUTTONS.ct2 = out.handles.s2;
    % add spheres if constrained inverse solution
    if isfield(model,'dipfit')...
            || ~isequal(model.xyz,zeros(1,3))
        try
            xyz = model.dipfit.Lpos;
            radius = model.dipfit.radius;
        catch
            xyz = model.xyz';
            radius = model.rad(1);
        end
        Np  = size(xyz,2);
        [x,y,z] = sphere(20);
        axes(D.PSD.handles.axes)
        for i=1:Np
            fvc = surf2patch(x.*radius+xyz(1,i),...
                y.*radius+xyz(2,i),z.*radius+xyz(3,i));
            D.PSD.handles.dipSpheres(i) = patch(fvc,...
                'parent',D.PSD.handles.axes,...
                'facecolor',[1 1 1],...
                'edgecolor','none',...
                'facealpha',0.5,...
                'tag','dipSpheres');
        end
        axis(D.PSD.handles.axes,'tight');
    end

    % plot time courses
    switch D.PSD.source.VIZU.timeCourses
        case 1
            Jp(1,:) = min(J,[],1);
            Jp(2,:) = max(J,[],1);
            D.PSD.source.VIZU.plotTC = plot(D.PSD.handles.axes2,...
                model.pst,Jp',...
                'color',0.5*[1 1 1],...
                'visible','off');
            % Add virtual electrode
            try
                ve = D.PSD.source.VIZU.ve;
            catch
                [mj ve] = max(max(abs(J),[],2));
                D.PSD.source.VIZU.ve =ve;
            end
            Jve = J(D.PSD.source.VIZU.ve,:);
            try
                qC  = model.qC(ve).*diag(model.qV)';
                ci  = 1.64*sqrt(qC);
                D.PSD.source.VIZU.pve2 = plot(D.PSD.handles.axes2,...
                    model.pst,Jve +ci,'b:',model.pst,Jve -ci,'b:');
            end
            D.PSD.source.VIZU.pve = plot(D.PSD.handles.axes2,...
                model.pst,Jve,...
                'color','b',...
                'visible','off');
        otherwise
            % this is meant to be extended for displaying something
            % else than just J (e.g. J^2, etc...)

    end
    D.PSD.source.VIZU.lineTime = line('parent',D.PSD.handles.axes2,...
        'xdata',[gridTime;gridTime],...
        'ydata',[miJ;maJ],...
        'visible','off');
    set(D.PSD.handles.axes2,...
        'ylim',[miJ;maJ]);  

    % create buttons
    object.type = 'buttons';
    object.list = [7;8;10];
    object.options.multSelect = 0;
    object.options.pst = pst;
    object.options.gridTime = gridTime;
    D = spm_eeg_review_uis(D,object);

    % create info text
    object.type = 'text';
    object.what = 'source';
    D = spm_eeg_review_uis(D,object);
    
    % set graphical object visible
    set(D.PSD.handles.mesh,'visible','on')
    set(D.PSD.handles.colorbar,'visible','on')
    set(D.PSD.handles.axes2,'visible','on')
    set(D.PSD.source.VIZU.lineTime,'visible','on')
    set(D.PSD.source.VIZU.plotTC,'visible','on')
    set(D.PSD.source.VIZU.pve,'visible','on')
    try
        set(D.PSD.handles.BMCplot,'visible','on');
        set(D.PSD.handles.hbar,'visible','on');
        set(D.PSD.handles.BMCcurrent,'visible','on');
        set(D.PSD.handles.BMCpanel,'visible','on');
    end
    

    set(D.PSD.handles.hfig,'userdata',D)



else

    uicontrol('style','text',...
        'units','normalized','Position',[0.14 0.84 0.7 0.04],...
        'string','There is no (imaging) inverse source reconstruction in this data file !',...
        'BackgroundColor',0.95*[1 1 1],...
        'tag','plotEEG')
    labels{1} = '1';
    callbacks{1} = [];
    hInv = D.PSD.handles.tabs.hp;
    spm_uitab(hInv,labels,callbacks,'plotEEG');


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GET DATA INFO
function [D] = DataInfo(D)


switch D.PSD.VIZU.uitable

    case 'off'

        % delete graphical objects from other main tabs
        delete(findobj('tag','plotEEG'));
        % create info text
        object.type = 'text';
        object.what = 'data';
        D = spm_eeg_review_uis(D,object);
        
        
        % add buttons
        object.type = 'buttons';
        object.list = [14,15];
        D = spm_eeg_review_uis(D,object);
        set(D.PSD.handles.BUTTONS.showSensors,...
            'position',[0.7 0.9 0.25 0.02]);
        set(D.PSD.handles.BUTTONS.saveHistory,...
            'string','save history as script',...
            'position',[0.7 0.87 0.25 0.02]);

    case 'on'



        if isempty(D.PSD.VIZU.fromTab) || ~isequal(D.PSD.VIZU.fromTab,'info')
            % delete graphical objects from other main tabs
            delete(findobj('tag','plotEEG'));
            % create info text
            object.type = 'text';
            object.what = 'data';
            D = spm_eeg_review_uis(D,object);
            % Create uitabs for channels and trials
            try
                D.PSD.VIZU.info;
            catch
                D.PSD.VIZU.info = 4;
            end
            labels = {'channels','trials','inv','history'};
            callbacks = {'spm_eeg_review_callbacks(''visu'',''main'',''info'',1)',...,...
                'spm_eeg_review_callbacks(''visu'',''main'',''info'',2)'...
                'spm_eeg_review_callbacks(''visu'',''main'',''info'',3)',...
                'spm_eeg_review_callbacks(''visu'',''main'',''info'',4)'};
            [h] = spm_uitab(D.PSD.handles.tabs.hp,labels,callbacks,'plotEEG',D.PSD.VIZU.info,0.9);
            D.PSD.handles.infoTabs = h;

        else
            % delete info table (if any)
            try delete(D.PSD.handles.infoUItable);end
            % delete info message (if any)
            try delete(D.PSD.handles.message);end
            % delete buttons if any
            try delete(D.PSD.handles.BUTTONS.OKinfo);end
            try delete(D.PSD.handles.BUTTONS.showSensors);end
            try delete(D.PSD.handles.BUTTONS.saveHistory);end
        end

        % add table and buttons
        object.type = 'buttons';
        object.list = [];

        switch D.PSD.VIZU.info

            case 1 % channels info
                object.list = [object.list;12;14];
                nc = length(D.channels);
                table = cell(nc,5);
                for i=1:nc
                    table{i,1} = D.channels(i).label;
                    table{i,2} = D.channels(i).type;
                    if D.channels(i).bad
                        table{i,3} = 'yes';
                    else
                        table{i,3} = 'no';
                    end
                    if ~isempty(D.channels(i).X_plot2D)
                        table{i,4} = 'yes';
                    else
                        table{i,4} = 'no';
                    end
                    table{i,5} = D.channels(i).units;
                end
                colnames = {'label','type','bad','position','units'};
                [ht,hc] = spm_uitable(table,colnames);
                set(ht,'units','normalized');
                set(hc,'position',[0.1 0.05 0.55 0.7],...
                    'tag','plotEEG');
                D.PSD.handles.infoUItable = ht;
                D.PSD.handles.infoUItable2 = hc;
                D = spm_eeg_review_uis(D,object); % this adds the buttons

            case 2 % trials info
                object.list = [object.list;12];
                ok = 1;
                if strcmp(D.type,'continuous')
                    try
                        ne = length(D.trials(1).events);
                        if ne == 0
                            ok = 0;
                        end
                    catch
                        ne = 0;
                        ok = 0;
                    end
                    if ne > 0
                        table = cell(ne,3);
                        for i=1:ne
                            table{i,1} = D.trials(1).label;
                            table{i,2} = D.trials(1).events(i).type;
                            table{i,3} = num2str(D.trials(1).events(i).value);
                            if ~isempty(D.trials(1).events(i).duration)
                                table{i,4} = num2str(D.trials(1).events(i).duration);
                            else
                                table{i,4} = [];
                            end
                            table{i,5} = num2str(D.trials(1).events(i).time);
                            table{i,6} = 'Undefined';
                            table{i,7} = num2str(D.trials(1).onset);
                        end
                        colnames = {'label','type','value','duration','time','bad','onset'};
                        [ht,hc] = spm_uitable(table,colnames);
                        set(ht,'units','normalized');
                        set(hc,'position',[0.1 0.05 0.74 0.7],...
                            'tag','plotEEG');
                    else
                        POS = get(D.PSD.handles.infoTabs.hp,'position');
                        D.PSD.handles.message = uicontrol('style','text','units','normalized',...
                            'Position',[0.14 0.84 0.7 0.04].*repmat(POS(3:4),1,2),...
                            'string','There is no event in this data file !',...
                            'BackgroundColor',0.95*[1 1 1],...
                            'tag','plotEEG');
                    end
                else
                    nt = length(D.trials);
                    table = cell(nt,3);
                    if strcmp(D.type,'single')
                        for i=1:nt
                            table{i,1} = D.trials(i).label;
                            ne = length(D.trials(i).events);
                            if ne == 0 || ((ne == 1) && isequal(D.trials(i).events(1).type, 'no events'))
                                table{i,2} = 'no events';
                                table{i,3} = 'no events';
                                table{i,4} = 'no events';
                                table{i,5} = 'no events';
                            elseif ne >1
                                table{i,2} = 'multiple events';
                                table{i,3} = 'multiple events';
                                table{i,4} = 'multiple events';
                                table{i,5} = 'multiple events';
                            else
                                table{i,2} = D.trials(i).events.type;
                                table{i,3} = num2str(D.trials(i).events.value);
                                if ~isempty(D.trials(i).events.duration)
                                    table{i,4} = num2str(D.trials(i).events.duration);
                                else
                                    table{i,4} = 'Undefined';
                                end
                                table{i,5} = num2str(D.trials(i).events.time);
                            end
                            if D.trials(i).bad
                                table{i,6} = 'yes';
                            else
                                table{i,6} = 'no';
                            end
                            table{i,7} = num2str(D.trials(i).onset);
                        end
                        colnames = {'label','type','value','duration','time','bad','onset'};
                        [ht,hc] = spm_uitable(table,colnames);
                        set(ht,'units','normalized');
                        set(hc,'position',[0.1 0.05 0.74 0.7],...
                            'tag','plotEEG');
                    else
                        for i=1:nt
                            table{i,1} = D.trials(i).label;
                            table{i,2} = num2str(D.trials(i).repl);
                            if D.trials(i).bad
                                table{i,3} = 'yes';
                            else
                                table{i,3} = 'no';
                            end
                        end
                        colnames = {'label','nb of repl','bad'};
                        [ht,hc] = spm_uitable(table,colnames);
                        set(ht,'units','normalized');
                        set(hc,'position',[0.1 0.05 0.32 0.7],...
                            'tag','plotEEG');
                    end
                end
                if ok
                    D.PSD.handles.infoUItable = ht;
                    D.PSD.handles.infoUItable2 = hc;
                    D = spm_eeg_review_uis(D,object); % this adds the buttons
                end

            case 3 % inv info

                object.list = [object.list;12];
                isInv = D.PSD.source.VIZU.isInv;
%                 isInv = 1:length(D.other.inv);
                if numel(isInv) >= 1 %D.PSD.source.VIZU.current ~= 0
                    Ninv = length(isInv);
                    table = cell(Ninv,12);
                    for i=1:Ninv
                        try
                            table{i,1} = [D.other.inv{isInv(i)}.comment{1},' '];
                        catch
                            table{i,1} = ' ';
                        end
                        table{i,2} = [D.other.inv{isInv(i)}.date(1,:)];
                        try
                            table{i,3} = [D.other.inv{isInv(i)}.inverse.modality];
                        catch
                            try
                                table{i,3} = [D.other.inv{isInv(i)}.modality];
                            catch
                                table{i,3} = '?';
                            end
                        end
                        table{i,4} = [D.other.inv{isInv(i)}.method];
                        try
                            table{i,5} = [num2str(length(D.other.inv{isInv(i)}.inverse.Is))];
                        catch
                            try
                                table{i,5} = [num2str(D.other.inv{isInv(i)}.inverse.n_dip)];
                            catch
                                table{i,5} = '?';
                            end
                        end
                        try
                            table{i,6} = [D.other.inv{isInv(i)}.inverse.type];
                        catch
                            table{i,6} = '?';
                        end
                        try
                            table{i,7} = [num2str(floor(D.other.inv{isInv(i)}.inverse.woi(1))),...
                                ' to ',num2str(floor(D.other.inv{isInv(i)}.inverse.woi(2))),' ms'];
                        catch
                            table{i,7} = [num2str(floor(D.other.inv{isInv(i)}.inverse.pst(1))),...
                                ' to ',num2str(floor(D.other.inv{isInv(i)}.inverse.pst(end))),' ms'];
                        end
                        try
                            if D.other.inv{isInv(i)}.inverse.Han
                                han = 'yes';
                            else
                                han = 'no';
                            end
                            table{i,8} = [han];
                        catch
                            table{i,8} = ['?'];
                        end
                        try
                            table{i,9} = [num2str(D.other.inv{isInv(i)}.inverse.lpf),...
                                ' to ',num2str(D.other.inv{isInv(i)}.inverse.hpf), 'Hz'];
                        catch
                            table{i,9} = ['?'];
                        end
                        try
                            table{i,10} = [num2str(size(D.other.inv{isInv(i)}.inverse.T,2))];
                        catch
                            table{i,10} = '?';
                        end
                        try
                            table{i,11} = [num2str(D.other.inv{isInv(i)}.inverse.R2)];
                        catch
                            table{i,11} = '?';
                        end
                        table{i,12} = [num2str(sum(D.other.inv{isInv(i)}.inverse.F))];
                    end
                    colnames = {'label','date','modality','model','#dipoles','method',...
                        'pst','hanning','band pass','#modes','%var','log[p(y|m)]'};
                    [ht,hc] = spm_uitable('set',table,colnames);
                    set(ht,'units','normalized');
                    set(hc,'position',[0.1 0.05 0.8 0.7],...
                        'tag','plotEEG');
                    D.PSD.handles.infoUItable = ht;
                    D.PSD.handles.infoUItable2 = hc;
                    D = spm_eeg_review_uis(D,object); % this adds the buttons
                else
                    POS = get(D.PSD.handles.infoTabs.hp,'position');
                    D.PSD.handles.message = uicontrol('style','text','units','normalized',...
                        'Position',[0.14 0.84 0.7 0.04].*repmat(POS(3:4),1,2),...
                        'string','There is no source reconstruction in this data file !',...
                        'BackgroundColor',0.95*[1 1 1],...
                        'tag','plotEEG');

                end


            case 4 % history info
                object.list = [object.list;15];
                
                table = spm_eeg_history(D);
                if ~isempty(table)
                    colnames = {'Process','function called','input file','output file'};
                    [ht,hc] = spm_uitable(table,colnames);
                    set(ht,'units','normalized','editable',0);
                    set(hc,'position',[0.1 0.05 0.8 0.7],...
                        'tag','plotEEG');
                    D.PSD.handles.infoUItable = ht;
                    D.PSD.handles.infoUItable2 = hc;
                else
                    POS = get(D.PSD.handles.infoTabs.hp,'position');
                    D.PSD.handles.message = uicontrol('style','text','units','normalized',...
                        'Position',[0.14 0.84 0.7 0.04].*repmat(POS(3:4),1,2),...
                        'string','The history of this file is not available !',...
                        'BackgroundColor',0.95*[1 1 1],...
                        'tag','plotEEG');
                end
                D = spm_eeg_review_uis(D,object); % this adds the buttons

        end

        % update data info if action called from 'info' tab...
        if ~isempty(D.PSD.VIZU.fromTab) && isequal(D.PSD.VIZU.fromTab,'info')
            [str] = spm_eeg_review_callbacks('get','dataInfo');
            set(D.PSD.handles.infoText,'string',str)
        end

end
