function DCM = spm_dcm_fnirs_specify(SPMf)
% Specify inputs of a DCM for fNIRS
% FORMAT [DCM] = spm_dcm_nirs_specify(SPMf)
%
% SPMf - SPM filename(s)
%
% DCM  - DCM structure (see spm_dcm_ui)
%__________________________________________________________________________
% Copyright (C) 2015-2016 Wellcome Trust Centre for Neuroimaging

% Will Penny & Sungho Tak
% $Id: spm_dcm_fnirs_specify.m 7060 2017-04-18 16:44:46Z will $

%--------------------------------------------------------------------------
%-Interactive window
Finter = spm_figure('GetWin','Interactive');
bcolor = get(Finter,'color');
WS     = spm('WinScale');
dx     = 20;

spm_input('Specify DCM for fNIRS:...  ',1,'d');

%--------------------------------------------------------------------------
% Get SPM.mat files
if ~nargin || isempty(SPMf)
    [SPMf, sts] = spm_select([1 Inf],'^SPM\.mat$','Select SPM.mat');
    if ~sts, DCM = []; return; end
end

SPMf = cellstr(SPMf);
n = numel(SPMf);
for i = 1:n
    load(SPMf{i});
    SPMn(i) = SPM;
end
clear SPM;

%--------------------------------------------------------------------------
%-Specify a name of DCM file and a directory to save DCM.mat
name  = spm_input('name for DCM_???.mat','+1','s');
[swd, sts] = spm_select(1, 'dir', 'Select a directory in which to place DCM-fNIRS results');
if ~sts, DCM = []; return; end

%--------------------------------------------------------------------------
%-Read optical measurements and their parameters
P.ns = [];
for i = 1:n
    fname     = SPMn(i).xY.VY;
    try
        raw   = load(fname, 'P');
    catch
        pth   = spm_file(SPMf{i},'fpath');
        fname = fullfile(pth, 'NIRS.mat');
        raw   = load(fname, 'P');
        raw.P.fname.pos = fullfile(pth, 'POS.mat');
    end
    P.ns(i)   = raw.P.ns;
    P.fname.nirs{i,1} = fname;
end
field = {'acoef', 'nch', 'wav', 'fs'};
for i = 1:size(field,2), P.(field{i}) = raw.P.(field{i}); end

%--------------------------------------------------------------------------
%-Read optical sensor positions
if isfield(raw.P.fname, 'pos')
    P.fname.pos = raw.P.fname.pos;
    clear raw;
else
    [P.fname.pos, sts] = spm_select(1, '^POS*.*\.mat$', 'Select an optical probe position file - POS.mat');
    if ~sts, DCM = []; return; end
end

%--------------------------------------------------------------------------
%-Preprocess time series of optical density changes (if necessary)
spm_input('Preprocess time series of optical density changes:...', 1, 'd');

% Temporal filtering
K.type = spm_input('temporal filtering?', '+1', 'Butterworth IIR|No');
if strcmpi(K.type, 'Butterworth IIR')
    K.cutoff = spm_input('stopband frequencies Hz [start end]:', '+1', 'r', '[0 0.008; 0.12 0.35; 0.7 1.5]');
    % defaults: [0 0.008]: very low-frequency drifts; [0.12 0.35]: respiration; [0.7 1.5]: cardiac pulsation
end

%--------------------------------------------------------------------------
%-Average time series over trials
W.type = spm_input('average time series over trials?', '+1', 'y/n');
if strcmpi(W.type, 'y')
    
    % read onset times from SPM.mat
    onsets = {}; names = {}; t0 = 0;
    for i = 1:n,
        nc = size(SPMn(i).Sess.U, 2);
        for j = 1:nc
            str = ['include ' SPMn(i).Sess.U(j).name{1} '?']; 
            if spm_input(str, '+1', 'y/n', [1 0], 1) 
                onsets = [onsets t0+SPMn(i).Sess.U(j).ons];
                names = [names SPMn(i).Sess.U(j).name];
            end
        end
        t0 = t0 + P.ns(i)./P.fs;
    end
    W.onsets = onsets;
    W.names = names;
    clear SPMn;
    
    spm_input('specify averaging window length [sec]:', '+1', 'd');
    nc = size(names, 2);
    fprintf('--------------------------------------------------------- \n');
    fprintf('time window [begin end]\n');
    fprintf('--------------------------------------------------------- \n');
    t0 = 0; wdur = min(diff(sort(cell2mat(W.onsets(:))))); 
    for i = 1:nc,
        W.durations(i) = spm_input(['for ' W.names{i} ':'], '+1', 'r', wdur, 1);
        fprintf('%-20s: %4.2f %4.2f [sec] \n', W.names{i}, t0, t0+W.durations(i));
        t0 = t0 + W.durations(i);
    end
    ns = sum(round(W.durations.*P.fs));
else
    ns = sum(P.ns);
end

%--------------------------------------------------------------------------
%-Specify regressors (eg, systemic confounds)
spm_input('Specify confounding effects:...', 1, 'd');
C.type = spm_input('regressors?', '+1', 'DCT|User|None');
switch C.type
    case 'DCT'
        C.period = spm_input('periods [sec]:', '+1','r', '[8 12]');
    case 'User'
        C.X0 = spm_input('enter matrix name or values:', '+1', 'r');
end

% Correction of pial vein contamination from fNIRS measurements
options.pialv = spm_input('correction of pial vein contamination?', '+1', 'b', {'yes', 'no'}, [1 0], 1);

%--------------------------------------------------------------------------
%-Specify fNIRS channels of interest
spm_input('Specify fNIRS channels to be used in DCM:...', 1, 'd');

% display channel positions
load(P.fname.pos);
spm_fnirs_viewer_sensor(R);

ans = spm_input('use all sensor measurements ?', '+1',  'y/n');
if strcmpi(ans, 'y')
    P.rois = R.ch.label; 
elseif strcmpi(ans, 'n')
    P.rois = spm_input('sensors of interest:', '+1', 'n');
end
clear R; 

%--------------------------------------------------------------------------
%-Specify hemodynamic source positions
spm_input('Specify hemo/neurodynamic source positions:...',1,'d');

m = spm_input('number of sources of interest', '+1', 'n');
xY = [];
for i = 1:m
    xY(i).name = spm_input(['name of source ' num2str(i)], '+1', 's');
    xY(i).xyz = spm_input('MNI coordinate [mm]', '+1', 'r')';
end

%--------------------------------------------------------------------------
%-Determine whether sources are spatially distributed
s = spm_input('type of hemodynamic source: ', '+1', 'b', {'Point','Distributed'}, [0 1], 2);
if ~s
    options.rs = 0;
else
    % maximum distance between point and distributed sources
    options.rs = spm_input('radius of distributed source [mm]', '+1', 'r', '4');
end

%--------------------------------------------------------------------------
%-Specify Green's function
[P.fname.g, sts] = spm_select(1, 'mat', 'Select a Greens function file');
if ~sts, DCM = []; return; end

%--------------------------------------------------------------------------
%-Store fNIRS variables in structure Y
Y.P = P; % parameters for measurements
Y.K = K; % temporal filtering
Y.W = W; % averaging window
Y.C = C; % confound regressor

%--------------------------------------------------------------------------
%-Specify experimental inputs
spm_input('Input specification:...  ',1,'d');

% specify onsets and durations of stimulus
UNITS = 'secs';
U   = {};
v   = spm_input('number of conditions/trials',2,'w1');

% get trials
for i = 1:v
    % get names
    str       = sprintf('name for condition/trial %d ?',i);
    Uname     = {spm_input(str,3,'s',sprintf('trial %d',i))};
    U(i).name = Uname;
    
    % get onsets
    str      = ['onsets - ' Uname{1} ' [sec]'];
    ons      = spm_input(str,4,'r',' ',[Inf 1]);
    U(i).ons = ons(:);
    
    % get durations
    str = 'duration[s] [sec]';
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
SPM.nscan = ns;
SPM.xBF.T = spm_get_defaults('stats.fmri.t'); % microtime resolution
SPM.xBF.dt = 1./(SPM.xBF.T * P.fs);
SPM.xBF.UNITS = 'secs';
SPM.Sess.U = U;
[Ut] = spm_get_ons(SPM,1);

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

%--------------------------------------------------------------------------
%- Model options (DCM-fNIRS uses default options)
options.two_state = 0;
options.nonlinear  = 0;
options.stochastic = 0;
options.centre     = 0;

%--------------------------------------------------------------------------
%- Graph connections
a     = zeros(m,m);
b     = zeros(m,m,nc);
c     = zeros(m,nc);
d     = zeros(m,m,0);

%--------------------------------------------------------------------------
%- Endogenous connections (A matrix)

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

%--------------------------------------------------------------------------
%-Get a
for i = 1:m
    for j = 1:m
        a(i,j) = get(h3(i,j),'Value');
    end
end

delete(findobj(get(Finter,'Children'),'flat'));

%--------------------------------------------------------------------------
%-Effects of causes (B and C matrices)
uicontrol(Finter,'String','done','Position', [300 100 060 020].*WS,...
    'Callback', 'uiresume(gcbf)');
for k = 1:nc
    
    %----------------------------------------------------------------------
    %-Buttons and labels
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
                
                %----------------------------------------------------------
                % Allow modulation of endogenous connections
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
    
    %----------------------------------------------------------------------
    %-Get c
    for i = 1:m
        c(i,k)   = get(h2(i),'Value');
    end
    
    %----------------------------------------------------------------------
    %-Get b allowing any 2nd order effects
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

%--------------------------------------------------------------------------
%-Effects of nonlinear modulations (D matrices)
if options.nonlinear
    uicontrol(Finter,'String','done','Position', [300 100 060 020].*WS,...
        'Callback', 'uiresume(gcbf)');
    for k = 1:m
        
        %------------------------------------------------------------------
        %-Buttons and labels
        str = sprintf('Effects of %-12s activity on connections',xY(k).name);
        spm_input(str,1,'d');
        
        for i = 1:m
            for j = 1:m
                if a(i,j)==1
                    
                    %------------------------------------------------------
                    % Allow modulation of endogenous connections
                    h4(i,j) = uicontrol(Finter,...
                        'Position',[220+dx*j 360-dx*i 020 020].*WS,...
                        'BackgroundColor',bcolor,...
                        'Style','radiobutton');
                end
            end
        end
        
        uiwait(Finter);
        
        %------------------------------------------------------------------
        %-Get d allowing any 2nd order effects
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

%--------------------------------------------------------------------------
%-Store all variables in DCM structure
DCM.Y = Y;
DCM.U = U;
DCM.a = a;
DCM.b  = b;
DCM.c = c;
DCM.d = d;
DCM.xY = xY;
DCM.v = ns;
DCM.n = length(xY); % number of sources of interest
DCM.options = options;

%-Save DCM.mat file
save(fullfile(swd, ['DCM_' name '.mat']), 'DCM', spm_get_defaults('mat.format'));
