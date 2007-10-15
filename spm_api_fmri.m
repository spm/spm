function varargout = spm_api_fmri(varargin)
% API for SPM.mat {fMRI}
% FORMAT spm_api_fmri(action,....)
%____________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_api_fmri.m 946 2007-10-15 16:36:06Z john $


% set action
%----------------------------------------------------------------------------
global study
try
    action = varargin{1};
catch
    action = 'initialise';
end

switch action
    
% open level 1 figure: spm_api_fmri('initialise')
%----------------------------------------------------------------------------
case {'initialise'}
    
    spm_defaults
    global defaults
    defaults.modality = 'FMRI';
    study = openfig('spm_study');
    global study
    h     = guihandles(study);
    set(h.swd,'String',pwd);

% delete subordinate figures: spm_api_fmri('delete',handle)
%----------------------------------------------------------------------------
case {'delete'}
    
    switch lower(get(varargin{2},'Tag'))
        
    case {'session','trial','parameter'}
       
        spm_api_fmri('quit',varargin{2});
        
    otherwise
 
        % remove from list and branch 
        %--------------------------------------------------------------------
        try
            h      = guihandles(varargin{2});
            i      = get(h.branch,'Value');
            P      = get(h.branch,'String');
            P(i)   = [];
            set(h.branch,'String',P)
            set(h.branch,'Value',1)
            branch = getappdata(h.branch,'branch');
            spm_api_fmri('delete',branch(i))
            branch(i) = [];
            setappdata(h.branch,'branch',branch)
        end
    end
    
% create a new figure: spm_api_fmri('new',handle,value)
%----------------------------------------------------------------------------
case {'new'}
    
    % get level
    %------------------------------------------------------------------------
    level = get(varargin{2},'Tag');
    h     = guihandles(varargin{2});
    
    switch level
        
    % new session
    %------------------------------------------------------------------------ 
    case {'study'}
        
        P = get(h.branch,'String');
        P{end + 1} = sprintf('session %i',length(P) + 1);
        set(h.branch,'String',P)
        fig    = openfig('spm_session');
        set(fig,'name',P{end});
        branch     = getappdata(h.branch,'branch');
        branch(end + 1) = fig;
        setappdata(h.branch,'branch',branch);
        
    % new trial
    %------------------------------------------------------------------------ 
    case {'session'}
        
        P = get(h.branch,'String');
        P{end + 1} = sprintf('trial %i',length(P) + 1);
        set(h.branch,'String',P)
        fig        = openfig('spm_trial');
        set(fig,'name',P{end});
        branch     = getappdata(h.branch,'branch');
        branch(end + 1) = fig;
        setappdata(h.branch,'branch',branch);
        
    % new parameter
    %------------------------------------------------------------------------ 
    case {'trial'}
        
        P  = get(h.branch,'String');
        P{end + 1} = sprintf('parameter %i',length(P) + 1);
        set(h.branch,'String',P)
        fig        = openfig('spm_parameter');
        set(fig,'name',P{end});
        branch     = getappdata(h.branch,'branch');
        branch(end + 1) = fig;
        setappdata(h.branch,'branch',branch);
        hp         = guihandles(fig);
        set(hp.P,'Userdata',get(h.ons,'Userdata'));
        spm_api_fmri('P',hp.P)
        
    end
    
% set values of SPM into figure: spm_api_fmri('set',handle,value)
%----------------------------------------------------------------------------
case {'set'}
    
    level = get(varargin{2},'Tag');
    h     = guihandles(varargin{2});
    
    switch level
        
    % level 1 - set entire study
    %---------------------------------------------------------------------
    case {'study'}
        
        SPM  = varargin{3};
        try, set(h.name, 'String', SPM.name), end
        try, set(h.study,'Name',   SPM.name), end
        try, set(h.SPMid,'String', SPM.SPMid),end
        
        % set level 1 values
        %---------------------------------------------------------------------
        set(h.RT,     'String',  sprintf('%.2f',SPM.xY.RT))
        set(h.RT,     'UserData',SPM.xY.RT)
        set(h.iGXcalc,'Value',  ~strcmp(lower(SPM.xGX.iGXcalc),'none'))
        set(h.form,   'Value',  ~strcmp(lower(SPM.xVi.form),'none'))
        set(h.HParam, 'Value',   isfinite(SPM.xX.K(1).HParam))
        
        set(h.UNITS,   'Value',  strcmp(SPM.xBF.UNITS,'secs') + 1)
        set(h.bf_name, 'Value',  SPM.xBF.order)
        set(h.order,   'Value',  SPM.xBF.order)
        set(h.length,  'Value',  fix(SPM.xBF.length/4))
        set(h.Volterra,'Value',  SPM.xBF.Volterra - 1)
        
        % set level subordinate levels
        %---------------------------------------------------------------------
        for i = 1:length(SPM.Sess)
            try
                name{i}  = SPM.Sess(i).name;
            catch
                name{i}  = sprintf('session %i',i);
            end
            
            % open, name and set children
            %-----------------------------------------------------------------
            fig(i) = openfig('spm_session');
            set(fig(i),'Name',name{i})
            spm_api_fmri('set',fig(i),SPM.Sess(i),SPM.xY.P(SPM.Sess(i).row,:))
        end
        
        % save handles of children in branch
        %---------------------------------------------------------------------
        branch = getappdata(h.branch,'branch');
        branch = [branch fig];
        setappdata(h.branch,'branch',branch)
        P      = get(h.branch,'String');
        if ~length(P), P = {}; end
        P      = {P{:}; name{:}};
        set(h.branch,'String',P)
        
    case {'session'}
        
        % set level 2 data (filenames)
        %---------------------------------------------------------------------
        Sess  = varargin{3};
        P     = varargin{4};
        set(h.P,'String', P)
        set(h.n,'String', num2str(size(P,1)))
        
        % covariates
        %---------------------------------------------------------------------
        C     = [];
        Cname = {};
        for i = 1:length(Sess.C)
            C     = [C Sess.C(i).C];
            Cname = {Cname{:} Sess.C(i).name{:}};
        end
        
        set(h.C,    'UserData', C)
        set(h.Cname,'String', Cname)
        spm_api_fmri('C',h.P)
        
        % set subordinate levels - trials
        %---------------------------------------------------------------------
        for i = 1:length(Sess.U)
            try
                name(i)  = Sess.U(i).name;
            catch
                name{i}  = sprintf('trial %i',i);
            end

            % open, name and set children
            %-----------------------------------------------------------------
            branch(i) = openfig('spm_trial');
            set(branch(i),'Name',name{i})
            spm_api_fmri('set',branch(i),Sess.U(i))
        end
        
        % save handles of children in branch
        %---------------------------------------------------------------------
        setappdata(h.branch,'branch',branch)  
        set(h.branch,'String',name)
        
    case {'trial'}
        
        % set level 3 data
        %---------------------------------------------------------------------
        U     = varargin{3};
        set(h.ons,'UserData',U.ons);
        set(h.ons,'String',sprintf('%0.2f ',U.ons));
        set(h.dur,'UserData',U.dur);
        set(h.dur,'String',sprintf('%0.2f ',U.dur));
        spm_api_fmri('u',h.u)
        
        % subordinate levels - trials
        %---------------------------------------------------------------------
        for i = 1:length(U.P)
            
            try
                name{i}  = U.P(i).name;
            catch
                name{i}  = sprintf('parameter %i',i);
            end
            
            % open, name and set children
            %-----------------------------------------------------------------
            branch(i) = openfig('spm_parameter');
            set(branch(i),'Name',name{i})
            spm_api_fmri('set',branch(i),U.P(i))
        end
        
        % save handles of children in branch
        %---------------------------------------------------------------------
        setappdata(h.branch,'branch',branch)  
        set(h.branch,'String',name)
        
    case {'parameter'}
        
        % set level 4 data
        %---------------------------------------------------------------------
        P     = varargin{3};
        try
            set(h.P,'UserData',P.P)
        catch
            set(h.P,'UserData',[])
        end
        try
            set(h.h,'Value',P.h + 1)
        catch
            set(h.h,'Value',1)
        end
        spm_api_fmri('P',h.P)
        if strcmp(get(h.parameter,'Name'),'none'), set(h.parameter,'Visible','off'),end
        
    end
        
case {'C'} 
        
    % get covariates
    %---------------------------------------------------------------------
    h = guihandles(varargin{2});
    C = get(h.C,'UserData');
    axes(h.C);cla
        
    % image covariates
    %---------------------------------------------------------------------
    n   = size(C,2);
    if n > 1
        imagesc(C)
        colormap(gray)
    end
    if n == 1
        plot(C)
    end
    axis off
        
    % check names
    %---------------------------------------------------------------------
    Cname = get(h.Cname,'String');
    m     = length(Cname);     
    if m > n
        Cname = Cname(1:n);
    elseif m < n
        for i = (m + 1):n
            Cname{end + 1} = sprintf('cov %i',i);
        end
    end
    set(h.Cname,'String',Cname)
        
case {'P'} 
        
    % get covariates
    %---------------------------------------------------------------------
    h = guihandles(varargin{2});
    P = get(h.P,'UserData');
    P = P(:);
    set(h.P,'UserData',P);
    axes(h.P);cla

    % image covariates
    %---------------------------------------------------------------------
    n   = length(P);
    if n > 1
        plot(P)
    end
    axis off
    
case {'u'} 
    
    % get stimulus function
    %---------------------------------------------------------------------
    h     = guihandles(varargin{2});
    axes(h.u);cla
    ons   = fix(get(h.ons,'Userdata')*16);
    dur   = fix(get(h.dur,'Userdata')*16);
    ons   = ons(:);
    if ~length(dur), dur = 0; end
    n     = length(ons);
    ons   = ons - min(ons) + 1;
    
    if length(dur) < n
        dur = sparse(1:n,1,dur(1),n,1);
    end
    dur   = dur(1:n);
    m     = ons(end) + dur(end) + 16;
    u     = sparse(ons,1,1,m,1) - sparse(ons + dur + 1,1,1,m,1);
    u     = cumsum(u);
    plot(full(u))
    axis off
    
    % number of trails
    %------------------------------------------------------------------------
    str = sprintf('%i onsets',n);
    set(h.n,'String',str)
   
% set values of SPM into figure: spm_api_fmri('get')
%----------------------------------------------------------------------------
case {'get'}
    
    global study
    h         = guihandles(study);
    SPM.name  = get(h.study,'name');
    
    SPM.SPMid        = get(h.SPMid,'String');
    SPM.swd          = get(h.swd,'String');
    SPM.xY.RT        = get(h.RT,   'Userdata');
    SPM.xGX.iGXcalc  = {'none'};
    SPM.xVi.form     = 'none';
    SPM.xX.K.HParam  = Inf;
    if get(h.iGXcalc,'Value'), SPM.xGX.iGXcalc = {'Scaling'}; end
    if get(h.form,   'Value'), SPM.xVi.form    = 'AR(1)'; end
    if get(h.HParam, 'Value'), SPM.xX.K.HParam = 128; end
    str              = get(h.UNITS,      'String');
    SPM.xBF.UNITS    = str{get(h.UNITS,  'Value')};
    str              = get(h.bf_name,    'String');
    SPM.xBF.name     = str{get(h.bf_name,'Value')};
    SPM.xBF.Volterra = get(h.Volterra,   'Value') + 1;
    SPM.xBF.order    = get(h.order,      'Value');
    SPM.xBF.length   = get(h.length,     'Value')*4;
    
    % get handles of children in branch (sessions)
    %-------------------------------------------------------------------------
    hsession = getappdata(h.branch,'branch');
    for s = 1:length(hsession)
        
        h      = guihandles(hsession(s));
        SPM.Sess(s).name = get(h.session,'name');
        
        P      = get(h.P,'String');
        C      = get(h.C,'UserData');
        Cname  = get(h.Cname,'String');
        nscan  = length(P);
        
        try
            for i = 1:size(P,1)
                c{i,1} = P(i,:);
            end
            P = c;
        end
        
        if size(C,1) ~= length(P) & size(C,1)
            figure(hsession(s))
            warndlg('check number of filenames and covariates')
            return
        end
        
        
        % filemames
        %---------------------------------------------------------------------
        try
            SPM.xY.P           = cat(1,SPM.xY.P,P);
            SPM.nscan          = [SPM.nscan nscan];
        catch
            SPM.xY.P           = P;
            SPM.nscan          = nscan;
        end
        
        % covariates
        %---------------------------------------------------------------------
        SPM.Sess(s).C.C    = C;
        SPM.Sess(s).C.name = Cname;
        
        % get handles of children in branch (trials)
        %---------------------------------------------------------------------
        htrial = getappdata(h.branch,'branch');
        for t = 1:length(htrial)
            h = guihandles(htrial(t));
            SPM.Sess(s).U(t).name   = {get(h.trial,'name')};  
            ons    = get(h.ons,'Userdata');
            dur    = get(h.dur,'Userdata');
            
            if length(dur) ~= 1 & length(dur) ~= length(ons)
                figure(htrial(t))
                warndlg('please check durations')
                return
            end
            
            % set onsets and durations
            %-----------------------------------------------------------------
            SPM.Sess(s).U(t).ons    = ons(:);
            SPM.Sess(s).U(t).dur    = dur(:);
            
            % get handles of children in branch (parameters)
            %-----------------------------------------------------------------
            hparam = getappdata(h.branch,'branch');
            if length(hparam)
                for p  = 1:length(hparam)
                    h  = guihandles(hparam(p));
                    SPM.Sess(s).U(t).P(p).name = get(h.parameter,'name');
                    P  = get(h.P,'Userdata');
                
                    if length(P) & (length(P) ~= length(SPM.Sess(s).U(t).ons))
                        figure(hparam(p))
                        warndlg('please check parameter vector')
                        return
                    end
                
                    SPM.Sess(s).U(t).P(p).P = P;
                    SPM.Sess(s).U(t).P(p).h = get(h.h,'Value') - 1;
                end
            else
                SPM.Sess(s).U(t).P.name = 'none';
                SPM.Sess(s).U(t).P.P = [];
                SPM.Sess(s).U(t).P.h = 0;
            end
        end
    end
    
    % assign output argument
    %------------------------------------------------------------------------
    try
        SPM.xY.P = strvcat(SPM.xY.P);
    end
    varargout{1} = spm_fMRI_design(SPM);
    
    
% delete figure and children: spm_api_fmri('quit',handle)
%----------------------------------------------------------------------------
case {'quit'}
    
    f = varargin{2};
    h = guihandles(f);
    try
        b = getappdata(h.branch,'branch');
    catch
        b = [];
    end
    for i = 1:length(b)
        spm_api_fmri('quit',b(i))
    end
    delete(f)
    
% re-set names: spm_api_fmri('refresh',handle)
%----------------------------------------------------------------------------
case {'refresh'}
    
    try
        f = varargin{2};
    catch
        global study, f = study;
    end
    
    h     = guihandles(f);
    try
        b = getappdata(h.branch,'branch');
    catch
        return
    end
    name  = {};
    for i = 1:length(b)
        name{i} = get(b(i),'Name');
        spm_api_fmri('refresh',b(i))
    end
    set(h.branch,'String',name)
    
% load a variable: x = spm_api_fmri('load')
%----------------------------------------------------------------------------
case {'load'}
    
    f     = spm_select(1,'.*t$','select MAT or text file');
    try
        x = load(f,'-ascii');
    catch
        x = load(f,'-mat');
        s = fieldnames(x);
        x = getfield(x,s{1});
    end
    varargout{1} = x;
    
% evalue string and place in UserData
%----------------------------------------------------------------------------
case {'eval'}

    h      = varargin{2};
    str    = get(h,'String');
    try, x = eval(['[' str ']']);    end
    try, x = eval(['[' str{:} ']']); end
    try
        set(h,'UserData',x)
        set(h,'String',sprintf('%.2f ',x));
    end
    
end
