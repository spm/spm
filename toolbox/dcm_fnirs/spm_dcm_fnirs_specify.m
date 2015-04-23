function DCM = spm_dcm_fnirs_specify(Y)
% Specify inputs of a DCM for fNIRS
% FORMAT [DCM] = spm_dcm_nirs_specify(Y)
%
% Y - structure array of fNIRS data (see spm_fnirs_convert.m) 
%
% DCM  - DCM structure (see spm_dcm_ui)
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny & Sungho Tak
% $Id: spm_dcm_fnirs_specify.m 6422 2015-04-23 16:55:52Z spm $

%-Interactive window
%--------------------------------------------------------------------------
Finter = spm_figure('GetWin','Interactive');
bcolor = get(Finter,'color');
WS     = spm('WinScale');
dx     = 20;

spm_input('Specify DCM for fNIRS:...  ',1,'d');

%-Name (DCM file) and its directory 
%--------------------------------------------------------------------------
name  = spm_input('name for DCM_???.mat','+1','s');

[swd, sts] = spm_select(1, 'dir', 'Select a directory where DCM.mat file will be written.');
if ~sts, DCM = []; return; end

%-Specify fNIRS data 
%--------------------------------------------------------------------------
if ~nargin || isempty(Y) 
    try 
        Y = spm_fnirs_convert; 
    catch 
        error('Cannot read fNIRS data file.');
    end
end

nscan = size(Y.y, 1);
dt = Y.dt;

%-Specify fNIRS channels of interest 
%--------------------------------------------------------------------------
str = 'fNIRS channels to be used for DCM analysis'; 
gui_s = spm_input(str, 1, 'b', {'All','Selected'}, [0 1], 2);
if ~gui_s % all channels 
    Y.ch = (1:Y.nch)';
else
    Y.ch = spm_input('Channels of interest:', '+1', 'n')';
    Y.nch = size(Y.ch, 1);
end
Y.y = reshape(Y.y(:, Y.ch, :), [nscan Y.nch * Y.nwav]); 

%-Specify source positions and 
% determine whether sources are spatially distributed
%--------------------------------------------------------------------------
m = spm_input('Number of sources of interest', '+1', 'n');
xY = [];
for i = 1:m
    xY(i).name = spm_input(['Name for source #' num2str(i)], '+1', 's');
    xY(i).xyz = spm_input('MNI coordinate [mm]', '+1', 'r')';
end

s = spm_input('Source: ', '+1', 'b', {'Point','Distributed'}, [0 1], 2);
if ~s
    options.rs = 0; 
else 
    options.rs = spm_input('Radius [mm]', '+1', 'r'); 
end

%-Specify confounding effects 
%--------------------------------------------------------------------------
K.type = spm_input('Confounding effects?', 1, 'DCT|Regressors|None'); 
switch K.type 
    case 'DCT'
        cutoff = spm_input('cutoff periods [sec]:', '+1','r', '[8 12]');
        K.cutoff = cutoff;
        K.row = 1:nscan;
        K.RT = dt;
    case 'Regressors' 
        C = spm_input('enter matrix name or values:', '+1', 'r'); 
        K.regressor = spm_detrend(C); 
end
options.pialv = spm_input('Pial vein oxygenation correction?', '+1', 'b', {'yes', 'no'}, [1 0], 1); 

%- Specify experimental inputs 
%--------------------------------------------------------------------------
spm_input('Input specification:...  ',1,'d');
% specify onsets and durations of stimulus 
UNITS = spm_input('specify design in','1','scans|secs');

U   = {};
v   = spm_input('number of conditions/trials',2,'w1');

% get trials
for i = 1:v
    % get names
    str       = sprintf('name for condition/trial %d ?',i);
    Uname     = {spm_input(str,3,'s',sprintf('trial %d',i))};
    U(i).name = Uname;
    
    % get onsets
    str      = ['vector of onsets - ' Uname{1}];
    ons      = spm_input(str,4,'r',' ',[Inf 1]);
    U(i).ons = ons(:);
    
    % get durations
    str = 'duration[s]';
    while 1
        dur = spm_input(str,5,'r',' ',[Inf 1]);
        if length(dur) == 1
            dur    = dur*ones(size(ons));
        end
        if length(dur) == length(ons), break, end
        str = sprintf('enter a scalar or [%d] vector',...
            length(ons));
    end
    U(i).dur = dur;
    
    U(i).P.name = 'none'; 
end
P.nscan = nscan;
P.xBF.T = spm_get_defaults('stats.fmri.t'); % microtime resolution 
P.xBF.dt = dt./P.xBF.T; 
P.xBF.UNITS = UNITS; 
P.Sess.U = U; 

[Ut] = spm_get_ons(P,1);

% with stimuli
U = [];
U.dt   =Ut(1).dt;
u      = length(Ut);
U.name = {};
U.u    = [];
for  i = 1:u
    for j = 1:length(Ut(i).name)
        U.u             = [U.u Ut(i).u(33:end,j)];
        U.name{end + 1} = Ut(i).name{j};
    end
end
nc = size(U.u, 2);

%- Model options (DCM-fNIRS uses default options) 
%--------------------------------------------------------------------------
options.two_state = 0; 
options.nonlinear  = 0;
options.stochastic = 0;
options.centre     = 0;

%- Graph connections 
%--------------------------------------------------------------------------
a     = zeros(m,m);
b     = zeros(m,m,nc);
c     = zeros(m,nc);
d     = zeros(m,m,0);

%- Endogenous connections (A matrix)
%--------------------------------------------------------------------------

%-Buttons and labels
%--------------------------------------------------------------------------
spm_input('Specify endogenous (fixed) connections from',1,'d')
spm_input('to',3,'d')
for i = 1:m
    str    = sprintf('%s %i',xY(i).name,i);
    h1(i)  = uicontrol(Finter,'String',str,...
        'Style','text',...
        'FontSize',10,...
        'BackgroundColor',bcolor,...
        'HorizontalAlignment','right',...
        'Position',[080 350-dx*i 080 020].*WS);
    h2(i)  = uicontrol(Finter,'String',sprintf('%i',i),...
        'Style','text',...
        'FontSize',10,...
        'BackgroundColor',bcolor,...
        'Position',[180+dx*i 350 010 020].*WS);
end
for i = 1:m
    for j = 1:m
        h3(i,j) = uicontrol(Finter,...
            'Position',[180+dx*j 350-dx*i 020 020].*WS,...
            'BackgroundColor',bcolor,...
            'Style','radiobutton');
        if i == j
            set(h3(i,j),'Value',1,...
                'enable','off');
        else
            set(h3(i,j),'enable','on','TooltipString', ...
                sprintf('from %s to %s',xY(j).name,xY(i).name));
        end
        if nc && i~=j
            set(h3(i,j),'Value',0);
        else
            set(h3(i,j),'Value',1);
        end
    end
end
uicontrol(Finter,'String','done','Position', [300 100 060 020].*WS,...
    'Callback', 'uiresume(gcbf)');

uiwait(Finter);

%-Get a
%--------------------------------------------------------------------------
for i = 1:m
    for j = 1:m
        a(i,j) = get(h3(i,j),'Value');
    end
end

delete(findobj(get(Finter,'Children'),'flat'));

%-Effects of causes (B and C matrices) 
%--------------------------------------------------------------------------
uicontrol(Finter,'String','done','Position', [300 100 060 020].*WS,...
    'Callback', 'uiresume(gcbf)');
for k = 1:nc
    
    %-Buttons and labels
    %----------------------------------------------------------------------
    str   = sprintf(...
        'Effects of %-12s on regions... and connections',...
        U.name{k});
    spm_input(str,1,'d');
    
    for i = 1:m
        h1(i)  = uicontrol(Finter,'String',xY(i).name,...
            'Style','text',...
            'BackgroundColor',bcolor,...
            'FontSize',10,...
            'Position',[080 350-dx*i 080 020].*WS);
        h2(i)  = uicontrol(Finter,...
            'Position',[160 360-dx*i 020 020].*WS,...
            'BackgroundColor',bcolor,...
            'Style','radiobutton');
    end
    for i = 1:m
        for j = 1:m
            if a(i,j) == 1
                
                % Allow modulation of endogenous connections
                %----------------------------------------------------------
                h3(i,j) = uicontrol(Finter,...
                    'Position',[220+dx*j 360-dx*i 020 020].*WS,...
                    'BackgroundColor',bcolor,...
                    'Style','radiobutton');
                set(h3(i,j),'TooltipString', ...
                    sprintf('from %s to %s',xY(j).name,xY(i).name));
                
            end
        end
    end
    
    uiwait(Finter);
    
    %-Get c
    %----------------------------------------------------------------------
    for i = 1:m
        c(i,k)   = get(h2(i),'Value');
    end
    
    %-Get b allowing any 2nd order effects
    %----------------------------------------------------------------------
    for i = 1:m
        for j = 1:m
            if a(i,j)==1
                b(i,j,k) = get(h3(i,j),'Value');
            end
        end
    end
    delete([h1(:); h2(:); h3(a==1)])
    
end
delete(findobj(get(Finter,'Children'),'flat'));


%-Effects of nonlinear modulations (D matrices)
%--------------------------------------------------------------------------
if options.nonlinear
    uicontrol(Finter,'String','done','Position', [300 100 060 020].*WS,...
        'Callback', 'uiresume(gcbf)');
    for k = 1:m
        
        %-Buttons and labels
        %------------------------------------------------------------------
        str = sprintf('Effects of %-12s activity on connections',xY(k).name);
        spm_input(str,1,'d');
        
        for i = 1:m
            for j = 1:m
                if a(i,j)==1
                    
                    % Allow modulation of endogenous connections
                    %------------------------------------------------------
                    h4(i,j) = uicontrol(Finter,...
                        'Position',[220+dx*j 360-dx*i 020 020].*WS,...
                        'BackgroundColor',bcolor,...
                        'Style','radiobutton');
                end
            end
        end
        
        uiwait(Finter);
        
        %-Get d allowing any 2nd order effects
        %------------------------------------------------------------------
        for i = 1:m
            for j = 1:m
                if a(i,j)==1
                    d(i,j,k) = get(h4(i,j),'Value');
                end
            end
        end
        delete(h4(a==1))
        
    end
end

delete(findobj(get(Finter,'Children'),'flat'));
spm_input('Thank you',1,'d')

%-Store all variables in DCM structure
%--------------------------------------------------------------------------
DCM.Y       = Y; 
DCM.U      = U; 
DCM.a       = a;
DCM.b       = b;
DCM.c       = c;
DCM.d       = d;
DCM.xY      = xY;
DCM.v       = nscan;
DCM.n       = length(xY); % number of sources of interest 
DCM.options = options;
DCM.K = K; 

%-Save
save(fullfile(swd, ['DCM_' name '.mat']), 'DCM', spm_get_defaults('mat.format'));
