function [] = spm_dcm_create (u)
% Specify a DCM model without having to use an SPM.mat file
% FORMAT [] = spm_dcm_create (u)
%
% u         eg. u(:,1), u(:,2) contain onsets for inputs
%           These can then be used at the user interface
% 
% This function is very much like spm_dcm_ui('specify') 
% but inputs etc. are specified via the user interface
%
% This function is intended to be used to create 
% DCM networks with known connectivity
% parameters from which synthetic data can be generated

Finter = spm_figure('GetWin','Interactive');
header = get(Finter,'Name');
set(Finter,'Name','Dynamic Causal Modelling')
WS     = spm('WinScale');

% name
%===================================================================
name  = spm_input('name for DCM_???.mat','+1','s');

% outputs
%===================================================================

% get cell array of region structures
%-------------------------------------------------------------------
% m, P, xY !!!

m = spm_input('Enter number of regions','+1','r',[],1);
for i=1:m,
    str=sprintf('Region %d',i);
    xY(i).name = spm_input(['Name for ',str],'+1','s');
    % Make up spurious VOI info
    % for compatability with spm_dcm_display
    xY(i).xyz = [i i i]'*10;
    xY(i).XYZmm = [i i i]'*10;
    xY(i).s=1;
    xY(i).spec=1;
end

% inputs
%===================================================================

% global parameters
global defaults
if ~isempty(defaults),
    fMRI_T  = defaults.stats.fmri.t;
    fMRI_T0 = defaults.stats.fmri.t0;
else,
    fMRI_T  = 16;
    fMRI_T0 = 1;
end;
SPM.xBF.T  = fMRI_T;
SPM.xBF.T0 = fMRI_T0;

spm_input('Basic parameters...',1,'d',mfilename)
SPM.xY.RT = spm_input('Interscan interval {secs}','+1','r',[],1);
SPM.nscan = spm_input(['scans per session e.g. 256'],'+1');
v=SPM.nscan;
SPM.xBF.dt = SPM.xY.RT/SPM.xBF.T;
str           = 'specify design in';
SPM.xBF.UNITS = spm_input(str,'+1','scans|secs');

Ui=spm_get_ons(SPM,1);

% Change input format to DCM input format
U.name = {};
U.u    = [];
for  i = 1:length(Ui)
    U.u             = [U.u Ui(i).u];
    U.name{end + 1} = Ui(i).name{1};
end
U.dt=Ui(1).dt;

% graph connections
%===================================================================
%n     = length(U);
n     = size(U.u,2);
a     = zeros(m,m);
c     = zeros(m,n);
b     = zeros(m,m,n);
d     = uicontrol(Finter,'String','done',...
    'Position',[300 50 060 020].*WS);
dx    = 35;
wx=30;
wy=20;

%WS = 2*WS;

%-intrinsic connections
%-------------------------------------------------------------------
spm_input('Specify intrinsic connections from',1,'d')
spm_input('to',3,'d')

for i = 1:m
    str    = sprintf('%s %i',xY(i).name,i);
    h1(i)  = uicontrol(Finter,'String',str,...
        'Style','text',...
        'HorizontalAlignment','right',...
        'Position',[080 356-dx*i 080 020].*WS);
    h2(i)  = uicontrol(Finter,'String',sprintf('%i',i),...
        'Style','text',...
        'Position',[180+dx*i 356 020 020].*WS);
end
for i = 1:m
    for j = 1:m
        cc=ceil([180+dx*j 360-dx*i wx wy].*WS);
        h3(i,j) = uicontrol(Finter,...
            'Position',cc,...
            'Style','edit');
        if i == j
            set(h3(i,j),'String','-1');
        else
            set(h3(i,j),'String','0');
        end
        
    end
end
drawnow

% wait for 'done'
%-----------------------------------------------------------
while(1)
    pause(0.01)
    if strcmp(get(gco,'Type'),'uicontrol')
        if strcmp(get(gco,'String'),'done')
            for i = 1:m
                for j = 1:m
                    a(i,j) = str2num(get(h3(i,j),'string'));
                end
            end
            delete([h1(:); h2(:); h3(:)])
            break
        end
    end
end


%-effects of causes
%-------------------------------------------------------------------
for k = 1:n
    
    % buttons and labels
    %-----------------------------------------------------------
%     str   = sprintf(...
%         'Effects of %-12s on regions... and connections',...
%         U.name{k});
    str   = sprintf(...
        'Effects of %-12s on regions... and connections',...
        U.name{k});
    spm_input(str,1,'d')
    
    
    for i = 1:m
        h1(i)  = uicontrol(Finter,'String',xY(i).name,...
            'Style','text',...
            'Position',[080 356-dx*i 080 020].*WS);
        h2(i)  = uicontrol(Finter,...
            'Position',[160 360-dx*i wx wy].*WS,...
            'Style','edit');
        set(h2(i),'String','0');
    end
    for i = 1:m
        for j = 1:m
            cc=ceil([200+dx*j 360-dx*i wx wy].*WS);
            h3(i,j) = uicontrol(Finter,...
                'Position',cc,...
                'Style','edit');
            set(h3(i,j),'String','0');
        end
    end
    drawnow
    
    % wait for 'done'
    %-----------------------------------------------------------
    set(gcf,'CurrentObject',h2(1))
    while(1)
        pause(0.01)
        if strcmp(get(gco,'Type'),'uicontrol')
            if strcmp(get(gco,'String'),'done')
                
                % get c
                %--------------------------------------------------
                for i = 1:m
                    c(i,k)   = str2num(get(h2(i),'string'));
                end
                
                % get b ensuring 2nd order effects are allowed
                %--------------------------------------------------
                for i = 1:m
                    for j = 1:m
                        b(i,j,k) = str2num(get(h3(i,j),'string'));
                        if i == j & ~c(i,k)
                            b(i,j,k) = 0;
                        end
                    end
                end
                delete([h1(:); h2(:); h3(:)])
                spm_input('Thank you',1,'d')
                break
                
            end
        end
    end
end
delete(d)

% Number of regions
n     = length(xY);


% ASSUME, for now, default hemodynamics

% Structure matrices
DCM.a=~(a==0);
DCM.b=~(b==0);
DCM.c=~(c==0);
[pE,pC,qE]=spm_dcm_priors(DCM.a,DCM.b,DCM.c);
DCM.H=qE;


% Now set up output structure 
%-------------------------------------------------------------------
X0    = ones(v,1);
Y.dt  = SPM.xY.RT;
Y.X0  = X0;
for i = 1:n,
    Y.name{i} = xY(i).name;
end
Y.Ce  = spm_Ce(ones(1,n)*v);

% Copy to data structure
DCM.A=a;
DCM.B=b;
DCM.C=c;
DCM.U=U;
DCM.Y=Y;
DCM.xY=xY;
DCM.v=v;
DCM.n=n;



%-Save and reset title
%-------------------------------------------------------------------
save(['DCM_',name],'DCM');

% Now generate output data
SNR=1;
spm_dcm_generate(DCM,SNR);

spm('FigName',header);
spm('Pointer','Arrow')

   
close 
