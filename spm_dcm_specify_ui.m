function DCM = spm_dcm_specify_ui(SPM,xY,settings)
% Interface for stepping the user through creating a DCM
% FORMAT DCM = spm_dcm_specify_ui(SPM,xY)
%
% SPM      - SPM structure from SPM.mat
% xY       - (optional) VOI structures to be inserted into the DCM
% settings - (optional) Structure of pre-populated settings for testing the
%            GUI without mouse clicks.
%
%             .u(i,j)      whether to include condition i regressor j
%             .delays      vector of delays [1 x n]
%             .TE          echo time
%             .nonlinear   non-linear DCM
%             .two_state   two-state DCM
%             .stochastic  stocastic DCM
%             .centre      mean-centring of inputs
%             .induced     induced responses)
%             .a .b .c .d  connectivity matrices
%
% DCM      - DCM structure (see spm_dcm_ui)
%__________________________________________________________________________
% Copyright (C) 2002-2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_specify_ui.m 7264 2018-02-22 14:43:47Z peter $


%-Interactive window
%--------------------------------------------------------------------------
f = struct();
f.Finter = spm_figure('GetWin','Interactive');
f.bcolor = get(f.Finter,'Color');
f.WS     = spm('WinScale');
f.dx     = 20;

spm_input('Specify DCM:...  ',1,'d');

%==========================================================================
% Outputs
%==========================================================================

%-Get cell array of region structures
%--------------------------------------------------------------------------
if nargin < 2 || isempty(xY)
    swd = SPM.swd;
    [P, sts] = spm_select([1 8],'^VOI.*\.mat$',{'select VOIs'},'',swd);
    if ~sts, DCM = []; return; end
    P = cellstr(P);
    
    xY = voi_files_to_array(P);
elseif iscell(xY) && ischar(xY{1})
    xY = voi_files_to_array(xY);  
end

m = numel(xY);

%==========================================================================
% Inputs
%==========================================================================

%-Get (nc) 'causes' or inputs U
%--------------------------------------------------------------------------
spm_input('Input specification:...  ',1,'d');
Sess   = SPM.Sess(xY(1).Sess);
if isempty(Sess.U)
    
    % spontaneous activity, i.e. no stimuli
    %----------------------------------------------------------------------
    U.u    = zeros(length(xY(1).u),1);
    U.name = {'null'};
    U.idx  = 0;
    
else
    
    % with stimuli
    %----------------------------------------------------------------------
    U.dt   = Sess.U(1).dt;
    u      = length(Sess.U);
    U.name = {};
    U.u    = [];
    U.idx  = [];

    for i = 1:u
        for j = 1:length(Sess.U(i).name)
            str = ['include ' Sess.U(i).name{j} '?'];
            
            try
                include_condition = (settings.u(i,j) == 1);
            catch
                include_condition = spm_input(str,'+1','y/n',[1 0],1);
            end
            
            if include_condition
                U.u             = [U.u Sess.U(i).u(33:end,j)];
                U.name{end + 1} = Sess.U(i).name{j};
                U.idx           = [U.idx; i j];
            end
        end
    end
    
    % Check for at least one (null) input
    %----------------------------------------------------------------------
    if isempty(U.u)
        U.u    = zeros(length(xY(1).u),1);
        U.name = {'null'};
        U.idx  = 0;
    end    
    
end

nc            = size(U.u,2);
is_endogenous = (nc == 1) && strcmp(U.name{1},'null');

%==========================================================================
% Timings
%==========================================================================

spm_input('Timing information:...  ',-1,'d');

%-VOI timings
%--------------------------------------------------------------------------
RT     = SPM.xY.RT;
t0     = spm_get_defaults('stats.fmri.t0');
t      = spm_get_defaults('stats.fmri.t');
T0     = RT * t0 / t;
try
    delays = settings.delays;
catch    
    delays = spm_input('VOI timings [s]','+1','r', repmat(T0,1,m),m,[0 RT]);
end

%-Echo time (TE) of data acquisition
%--------------------------------------------------------------------------
TE    = 0.04;
TE_ok = 0;

str = { 'Extreme value for TE or TE undefined.',...
                'Please re-enter TE (in seconds!)'};
    
try
    TE = settings.TE;    
catch        
    while ~TE_ok    
        TE = spm_input('Echo time, TE [s]', '+1', 'r', TE);

        if ~TE || (TE < 0) || (TE > 0.1)            
            spm_input(str,'+1','bd','OK',[1],1);
        else
            TE_ok = 1;
        end
    end
end

%==========================================================================
% Model options
%==========================================================================
spm_input('Model options:...  ',-1,'d');
try
    options.nonlinear = settings.nonlinear;
catch
    options.nonlinear = spm_input('modulatory effects','+1','b',{'bilinear','nonlinear'},[0 1],1);
end

try
    options.two_state = settings.two_state;
catch
    options.two_state = spm_input('states per region', '+1','b',{'one','two'},[0 1],1);
end

try
    options.stochastic = settings.stochastic;
catch
    options.stochastic = spm_input('stochastic effects','+1','b',{'no','yes'}, [0 1],1);
end

try
    options.centre = settings.centre;
catch
    options.centre = spm_input('centre input', '+1','b',{'no','yes'}, [0 1],1);
end

if options.stochastic 
    options.induced = 0;
else
    try
        options.induced = settings.induced;
    catch
        options.induced = spm_input('fit timeseries or CSD','+1','b',{'timeseries','CSD'}, [0 1],1);    
    end
end

%==========================================================================
% Graph connections
%==========================================================================

%-Endogenous connections (A matrix)
try
    a = settings.a;
catch
    a = query_user_for_a(m, xY, is_endogenous, f);
end

%-Effects of causes (B and C matrices)
try
    b = settings.b;
    c = settings.c;
catch
    [b,c] = query_user_for_bc(a, xY, U, nc, is_endogenous, f);
end

%-Effects of nonlinear modulations (D matrices)
try
    d = settings.d;
catch
    d = query_user_for_d(a, xY, options.nonlinear, f);
end

% Finish
spm_input('Thank you',1,'d')

%==========================================================================
% Prepare DCM
%==========================================================================

%-Response variables & confounds (NB: the data have been whitened)
%--------------------------------------------------------------------------
n     = length(xY);                      % number of regions
v     = length(xY(1).u);                 % number of time points
Y.dt  = SPM.xY.RT;
Y.X0  = xY(1).X0;
for i = 1:n
    Y.y(:,i)  = xY(i).u;
    Y.name{i} = xY(i).name;
end

%-Error precision components (one for each region) - i.i.d. (because of W)
%--------------------------------------------------------------------------
Y.Q        = spm_Ce(ones(1,n)*v);

%-Store all variables in DCM structure
%--------------------------------------------------------------------------
DCM.a       = a;
DCM.b       = b;
DCM.c       = c;
DCM.d       = d;
DCM.U       = U;
DCM.Y       = Y;
DCM.xY      = xY;
DCM.v       = v;
DCM.n       = n;
DCM.TE      = TE;
DCM.delays  = delays;
DCM.options = options;

%==========================================================================
function a = query_user_for_a(m,xY,is_endogenous,f)
% Query user for the A-matrix (endogenous connections)
%
% m             - number of regions
% xY            - regions of interest
% is_endogenous - if true, there are no inputs (resting state)
% f             - struct of figure variables

a = zeros(m,m);

%-Buttons and labels
spm_input('Specify endogenous (fixed) connections from',1,'d')
spm_input('to',3,'d')
for i = 1:m
    str    = sprintf('%s %i',xY(i).name,i);
    h1(i)  = uicontrol(f.Finter,'String',str,...
        'Style','text',...
        'FontSize',10,...
        'BackgroundColor',f.bcolor,...
        'HorizontalAlignment','right',...
        'Position',[080 350-f.dx*i 080 020].*f.WS);
    h2(i)  = uicontrol(f.Finter,'String',sprintf('%i',i),...
        'Style','text',...
        'FontSize',10,...
        'BackgroundColor',f.bcolor,...
        'Position',[180+f.dx*i 350 010 020].*f.WS);
end
for i = 1:m
    for j = 1:m
        h3(i,j) = uicontrol(f.Finter,...
            'Position',[180+f.dx*j 350-f.dx*i 020 020].*f.WS,...
            'BackgroundColor',f.bcolor,...
            'Style','radiobutton');
        if i == j
            set(h3(i,j),'Value',1,...
                'enable','off');
        else
            set(h3(i,j),'enable','on','TooltipString', ...
                sprintf('from %s to %s',xY(j).name,xY(i).name));
        end
        if ~is_endogenous && i~=j
            set(h3(i,j),'Value',0);
        else
            set(h3(i,j),'Value',1);
        end
    end
end
uicontrol(f.Finter,'String','done','Position', [300 100 060 020].*f.WS,...
    'Callback', 'uiresume(gcbf)');

uiwait(f.Finter);

%-Get a
for i = 1:m
    for j = 1:m
        a(i,j) = get(h3(i,j),'Value');
    end
end

delete(findobj(get(f.Finter,'Children'),'flat'));

uicontrol(f.Finter,'String','done','Position', [300 100 060 020].*f.WS,...
    'Callback', 'uiresume(gcbf)');


%==========================================================================
function [b,c] = query_user_for_bc(a, xY, U, nc, is_endogenous, f)
% Query user for the B- and C-matrices (network inputs)
%
% a             - connectivity matrix
% xY            - regions of interest
% U             - input structure
% nc            - number of causes (inputs)
% is_endogenous - if true, there are no inputs (resting state)
% f             - struct of figure variables

m = size(a,1);

b = zeros(m,m,nc);
c = zeros(m,nc);

if is_endogenous, return; end

for k = 1:nc

    %-Buttons and labels
    str   = sprintf(...
        'Effects of %-12s on regions... and connections',...
        U.name{k});
    spm_input(str,1,'d');

    for i = 1:m
        h1(i)  = uicontrol(f.Finter,'String',xY(i).name,...
            'Style','text',...
            'BackgroundColor',f.bcolor,...
            'FontSize',10,...
            'Position',[080 350-f.dx*i 080 020].*f.WS);
        h2(i)  = uicontrol(f.Finter,...
            'Position',[160 360-f.dx*i 020 020].*f.WS,...
            'BackgroundColor',f.bcolor,...
            'Style','radiobutton');
    end
    for i = 1:m
        for j = 1:m
            if a(i,j) == 1

                % Allow modulation of endogenous connections
                h3(i,j) = uicontrol(f.Finter,...
                    'Position',[220+f.dx*j 360-f.dx*i 020 020].*f.WS,...
                    'BackgroundColor',f.bcolor,...
                    'Style','radiobutton');
                set(h3(i,j),'TooltipString', ...
                    sprintf('from %s to %s',xY(j).name,xY(i).name));

            end
        end
    end

    uiwait(f.Finter);

    %-Get c    
    for i = 1:m
        c(i,k)   = get(h2(i),'Value');
    end

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

delete(findobj(get(f.Finter,'Children'),'flat'));

%==========================================================================
function d = query_user_for_d(a, xY, is_nonlinear, f)
% Query user for the D-matrix (non-linear inputs)
%
% a             - endogenous connectivity matrix
% xY            - regions of interest
% is_nonlinear  - if true, this is a non-linear DCM
% f             - struct of figure variables

m = size(a,1);
d = zeros(m,m,0);

if ~is_nonlinear, return; end

uicontrol(f.Finter,'String','done','Position', [300 100 060 020].*f.WS,...
    'Callback', 'uiresume(gcbf)');
for k = 1:m

    %-Buttons and labels
    str = sprintf('Effects of %-12s activity on connections',xY(k).name);
    spm_input(str,1,'d');

    for i = 1:m
        for j = 1:m
            if a(i,j)==1

                % Allow modulation of endogenous connections
                h4(i,j) = uicontrol(f.Finter,...
                    'Position',[220+f.dx*j 360-f.dx*i 020 020].*f.WS,...
                    'BackgroundColor',f.bcolor,...
                    'Style','radiobutton');
            end
        end
    end

    uiwait(f.Finter);

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

delete(findobj(get(f.Finter,'Children'),'flat'));

%==========================================================================
function xY = voi_files_to_array(P)
% Convert a cell array of VOI files to an array
m  = numel(P);
xY = [];
for i = 1:m
    p  = load(P{i},'xY');
    xY = spm_cat_struct(xY,p.xY);
end