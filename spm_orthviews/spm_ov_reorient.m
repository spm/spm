function ret = spm_ov_reorient(varargin)
% Reorient tool - plugin for spm_orthviews
%
% This tool provides the capabilities of the reorientation widget in SPM's
% "Display" for any image displayed within spm_orthviews.  The control
% fields are drawn in the SPM Interactive window and work as described in
% the Display routine.
% The advantage of using this tool within CheckReg is that it allows to
% reorient images while comparing their position to reference images
% simultaneously.
%
% This routine is a plugin to spm_orthviews. For general help about
% spm_orthviews and plugins type
%             help spm_orthviews
% at the MATLAB prompt.
%__________________________________________________________________________

% Volkmar Glauche
% Copyright (C) 2012-2022 Wellcome Centre for Human Neuroimaging


global st
if isempty(st)
    error('This routine can only be called as a plugin for spm_orthviews.');
end

if nargin < 2
    error(['Not enough input arguments.\n',...
           'Usage: spm_orthviews(''reorient'', cmd, volhandle, varargin)']);
end

cmd = lower(varargin{1});
volhandle = varargin{2};

switch cmd
    
    %-Context menu and callbacks
    %----------------------------------------------------------------------
    case 'context_menu'
        item0 = uimenu(varargin{3}, 'Label', 'Reorient',...
            'Tag', ['REORIENT_M_', num2str(volhandle)]);
        uimenu(item0, 'Label', 'This image', 'Callback', ...
            ['spm_orthviews(''reorient'',''context_init'', ', ...
            num2str(volhandle), ');'],...
            'Tag', ['REORIENT_0_', num2str(volhandle)]);
        uimenu(item0, 'Label', 'All images', 'Callback', ...
            ['spm_orthviews(''reorient'',''context_init'', 0);'],...
            'Tag', ['REORIENT_0_', num2str(volhandle)]);
        uimenu(item0, 'Label', 'Set origin to crosshair', 'Callback', ...
            ['spm_orthviews(''reorient'',''context_origin'', ', ...
            num2str(volhandle), ');'],...
            'Tag', ['REORIENT_0_', num2str(volhandle)]);
        uimenu(item0, 'Label', 'Quit Reorient image', ...
            'Tag', ['REORIENT_1_', num2str(volhandle)], ...
            'Visible', 'off');
        uimenu(item0, 'Label', 'Help', 'Callback', ...
            sprintf('spm_help(''%s'');', mfilename));
        ret = item0;
        
    case 'context_init'
        Finter = spm_figure('GetWin', 'Interactive');
        Fgraph = spm_figure('GetWin', 'Graphics');
        figure(Finter);
        spm_input('!DeleteInputObj',Finter);
        handles = spm_orthviews('valid_handles');
        labels = {'right  {mm}', 'forward  {mm}', 'up  {mm}',...
            'pitch  {rad}', 'roll  {rad}', 'yaw  {rad}',...
            'resize  {x}', 'resize  {y}', 'resize {z}'};
        tooltips = {'Translation', 'Translation', 'Translation', ...
            'Rotation', 'Rotation', 'Rotation', ...
            'Zoom', 'Zoom', 'Zoom', ...
            'Number of contour lines'};
        if volhandle == 0 || numel(handles) == 1
            % Reorient all images
            volhandle = handles;
            prms(7:9) = 1;
        else
            % get initial parameter values from st.vols{volhandle}.premul
            prms = spm_imatrix(st.vols{volhandle}.premul);
            prms(10) = 3; % default number of contour lines
            labels{end+1} = 'contour lines';
        end
        
        WS = spm('WinScale');
        u0 = uipanel(Finter, 'Units','Pixels', 'Title','',...
            'Position',[75 60 255 275].*WS);
        u1 = uipanel('Parent',u0, 'Units','Pixels', 'Title','',...
            'Position',[5 5 244 264].*WS,...
            'BorderType','Line', 'HighlightColor',[0 0 0]);
        st.vols{volhandle(1)}.reorient.p = [u0, u1];
        
        st.vols{volhandle(1)}.reorient.order = uicontrol('Parent',u1,...
            'Style','PopupMenu', 'Position', [5 32.5 232 025].*WS, ...
            'BackgroundColor',0.94*[1 1 1],...
            'String',{'Translation(1) Rotation(2) Zoom(3)', ...
            'Zoom(1) Translation(2) Rotation(3)', ...
            'Zoom(1) Rotation(2) Translation(3)'},...
            'Callback',['spm_orthviews(''reorient'',''reorient'',[',...
            num2str(volhandle),']);']);
        st.vols{volhandle(1)}.reorient.b(1) = uicontrol('Parent',u1,...
            'Style','PushButton', 'Position',[40 5 165 025].*WS, ...
            'BackgroundColor',0.94*[1 1 1], ...
            'String','Apply to image(s)', ...
            'Callback',['spm_orthviews(''reorient'',''apply'',[',...
            num2str(volhandle), ']);']);
        for k = handles
            % Find context menu handles
            obj = findobj(Fgraph, 'Tag',  ['REORIENT_M_', num2str(k)]);
            if any(k == volhandle)
                % Show 'Quit Reorient' for images being reoriented
                objh = findobj(obj, 'Tag', ['REORIENT_0_', num2str(k)]);
                objs = findobj(obj, 'Tag', ['REORIENT_1_', num2str(k)]);
                set(objh,'Visible','off');
                set(objs, 'Callback', ...
                    ['spm_orthviews(''reorient'',''context_quit'', [', ...
                    num2str(volhandle), ']);'],'Visible','on');
                st.vols{k}.reorient.oldpremul = st.vols{k}.premul;
            else
                % Do not show 'Reorient Images' context menu in other images
                set(obj, 'Visible', 'off');
            end
        end
        hpos = 240:-20:60;
        for k = 1:numel(labels)
            st.vols{volhandle(1)}.reorient.l(k) = uicontrol('Parent',u1,...
                'Style','Text', 'BackgroundColor',0.94*[1 1 1], ...
                'Position',[40 hpos(k) 100 016].*WS, 'String',labels{k});
            st.vols{volhandle(1)}.reorient.e(k) = uicontrol('Parent',u1,...
                'Style','edit', 'BackgroundColor',[1 1 1], ...
                'Callback',['spm_orthviews(''reorient'',''reorient'',[',...
                num2str(volhandle),'])'], ...
                'Position',[140 hpos(k) 065 020].*WS, 'String',num2str(prms(k)), ...
                'ToolTipString',tooltips{k});
        end
        spm_orthviews('redraw');
        
    case 'context_quit'
        Finter = spm_figure('FindWin', 'Interactive');
        Fgraph = spm_figure('FindWin', 'Graphics');
        try
            delete(st.vols{volhandle(1)}.reorient.e);
            delete(st.vols{volhandle(1)}.reorient.l);
            delete(st.vols{volhandle(1)}.reorient.b);
            delete(st.vols{volhandle(1)}.reorient.order);
            delete(st.vols{volhandle(1)}.reorient.p);
        end
        try
            for v = volhandle
                if ~isempty(st.vols{v}.reorient.lh)
                    delete(cat(1,st.vols{v}.reorient.lh{:}));
                end
            end
        end
        
        for k = spm_orthviews('valid_handles')
            try
                st.vols{k}.premul = st.vols{k}.reorient.oldpremul;
                st.vols{k} = rmfield(st.vols{k},'reorient');
            end
            obj = findobj(Fgraph, 'Tag',  ['REORIENT_M_', num2str(k)]);
            if any(k == volhandle)
                objh = findobj(obj, 'Tag', ['REORIENT_1_', num2str(k)]);
                set(objh, 'Visible', 'off', 'Callback','');
                objs = findobj(obj, 'Tag', ['REORIENT_0_', num2str(k)]);
                set(objs, 'Visible', 'on');
            else
                set(obj, 'Visible', 'on');
            end
        end
        spm_orthviews('redraw');
        
    case 'context_origin'
        pos = spm_orthviews('pos');
        M = spm_matrix(-pos');
        P = {spm_file(st.vols{volhandle}.fname, 'number', st.vols{volhandle}.n)};
        p = spm_fileparts(st.vols{volhandle}.fname);
        [P, sts] = spm_select(Inf, 'image', {'Image(s) to reorient'}, P, p);
        if ~sts
            disp('''Set origin to crosshair'' cancelled.');
            return;
        end
        sv = questdlg('Save reorientation matrix for future reference?',...
            'Save Matrix','Yes','No','Yes');
        if strcmpi(sv, 'yes')
            [p,n] = spm_fileparts(st.vols{volhandle}.fname);
            [f,p] = uiputfile(fullfile(p, [n '_reorient.mat']));
            if ~isequal(f,0)
                save(fullfile(p,f),'M', spm_get_defaults('mat.format'));
            end
        end
        if ~isempty(P)
            spm_jobman('serial', '', 'spm.util.reorient', cellstr(P), M);
            spm_orthviews('reload_mats');
            spm_orthviews('Reposition', [0 0 0]);
        end
        
    % Interaction callbacks
    %----------------------------------------------------------------------
        
    case 'apply'
        M = st.vols{volhandle(1)}.premul;
        N = numel(volhandle);
        P = cell(N, 1);
        for i = 1:N
            P{i} = spm_file(st.vols{volhandle(i)}.fname, ...
                'number', st.vols{volhandle(i)}.n);
        end
        [P, sts] = spm_select(Inf, 'image', {'Image(s) to reorient'}, P);
        if ~sts
            disp('Reorientation cancelled.');
            spm_orthviews('reorient', 'context_quit', volhandle);
            return;
        end
        sv = questdlg('Save reorientation matrix for future reference?', ...
            'Save Matrix','Yes','No','Yes');
        if strcmpi(sv, 'yes')
            [p,n] = spm_fileparts(st.vols{volhandle(1)}.fname);
            [f,p] = uiputfile(fullfile(p, [n '_reorient.mat']));
            if ~isequal(f,0)
                save(fullfile(p,f),'M', spm_get_defaults('mat.format'));
            end
        end        
        if isempty(P), return, end
        spm_jobman('serial', '', 'spm.util.reorient', cellstr(P), M);
        for i = 1:N
            st.vols{volhandle(i)}.reorient.oldpremul = eye(4);
        end
        spm_orthviews('reload_mats');
        spm_orthviews('reorient', 'context_quit', volhandle);
        
    case 'reorient'
        prms=zeros(1,12);
        for k=1:9
            prms(k) = evalin('base',[get(st.vols{volhandle(1)}.reorient.e(k),'string') ';']);
        end
        switch get(st.vols{volhandle(1)}.reorient.order, 'value')
            case 1
                order = 'Z*S*R*T';
            case 2
                order = 'R*T*Z*S';
            case 3
                order = 'T*R*Z*S';
        end
        for k = volhandle
            st.vols{k}.premul = spm_matrix(prms,order) * ...
                st.vols{k}.reorient.oldpremul;
        end
        spm_orthviews('redraw');
        
    case 'redraw'
        if isfield(st.vols{volhandle}.reorient, 'e') && numel(st.vols{volhandle}.reorient.e)==10
            if isfield(st.vols{volhandle}.reorient,'lh')
                if ~isempty(st.vols{volhandle}.reorient.lh)
                    delete(cat(1,st.vols{volhandle}.reorient.lh{:}));
                end
            end
            st.vols{volhandle}.reorient.lh = {};
            ncl = str2double(get(st.vols{volhandle}.reorient.e(10),'string'));
            if ncl > 0
                todraw = setxor(spm_orthviews('valid_handles'),volhandle);
                for d = 1:3
                    CData = sqrt(sum(get(st.vols{volhandle}.ax{d}.d,'CData').^2, 3));
                    CData(isinf(CData)) = NaN;
                    for h = todraw
                        set(st.vols{h}.ax{d}.ax,'NextPlot','add');
                        [C,st.vols{volhandle}.reorient.lh{end+1}] = ...
                            contour(st.vols{h}.ax{d}.ax,CData,ncl,'r-');
                    end
                end
                set(cat(1,st.vols{volhandle}.reorient.lh{:}),'HitTest','off');
            end
        end
        
    otherwise
        warning('Unknown action: %s', cmd);
end
