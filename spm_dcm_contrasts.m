function [con] = spm_dcm_contrasts (DCM_filename,D)
% Make contrast vector for a DCM model 
% FORMAT [con] = spm_dcm_contrasts (DCM_filename,D)
%
% DCM_filename  DCM file name
% D             'A','B' or 'C' ie. contrast for which connectivity matrix
%
% con           Column vector specifying contrast of parameters
%
% %W% Will Penny %E%

Finter = spm_figure('GetWin','Interactive');
header = get(Finter,'Name');
set(Finter,'Name','Dynamic Causal Modelling')
WS     = spm('WinScale');

P{1} = DCM_filename;
load(P{:})

m=DCM.n;
xY=DCM.xY;
U=DCM.U;

n     = size(DCM.U.u,2);
a     = zeros(m,m);
c     = zeros(m,n);
b     = zeros(m,m,n);

if (nargin < 2) | isempty(D)
    str     = 'contrast for';
    D       = spm_input(str,1,'b',{'A','B','C'});
end


dx    = 35;
wx=30;
wy=20;

d     = uicontrol(Finter,'String','done',...
    'Position',[300 50 060 020].*WS);

switch D,
    
case 'A' % intrinsic
    
    str   = sprintf(...
            'Enter contrast for A: ');
        spm_input(str,1,'d');
        
    % Print names and numbers of regions
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
    
    % Set contrast values to zero and display
    for i = 1:m
        for j = 1:m
            cc=ceil([180+dx*j 360-dx*i wx wy].*WS);
            h3(i,j) = uicontrol(Finter,...
                'Position',cc,...
                'Style','edit');
            set(h3(i,j),'String','0');
            
        end
    end
    drawnow
    
    % wait for use to i/p contrast and press 'done'
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
                spm_input(' ',1,'d')
                break
            end
        end
    end
    delete(d)
    con=a(:);
    
    
case 'B' % modulatory
    %---------------------------------------------------
    
    for k = 1:n,
        str   = sprintf(...
            'Enter contrast for B: effects of %-12s',...
            U.name{k});
        spm_input(str,1,'d')
        
        for i = 1:m
            h1(i)  = uicontrol(Finter,'String',xY(i).name,...
                'Style','text',...
                'Position',[080 356-dx*i 080 020].*WS);
        end
        for i = 1:m
            for j = 1:m
                cc=ceil([180+dx*j 360-dx*i wx wy].*WS);
                h3(i,j) = uicontrol(Finter,...
                    'Position',cc,...
                    'Style','edit');
                set(h3(i,j),'String','0');
            end
        end
        drawnow
        
        % wait for 'done'
        %-----------------------------------------------------------
        set(gcf,'CurrentObject',h3(1))
        while(1)
            pause(0.01)
            if strcmp(get(gco,'Type'),'uicontrol')
                if strcmp(get(gco,'String'),'done')
                    for i = 1:m
                        for j = 1:m
                            b(i,j,k) = str2num(get(h3(i,j),'string'));
                        end
                    end
                    delete([h1(:); h3(:)])
                    spm_input(' ',1,'d')
                    break
                    
                end
            end
        end
        
    end
    delete(d)
    con=b(:);
    
    
case 'C' % input
    %---------------------------------------------------
    
    for k = 1:n,
        str   = sprintf(...
            'Enter contrast for C: Effects of %-12s',...
            U.name{k});
        spm_input(str,1,'d');
        
        for i = 1:m
            h1(i)  = uicontrol(Finter,'String',xY(i).name,...
                'Style','text',...
                'Position',[080 356-dx*i 080 020].*WS);
            h2(i)  = uicontrol(Finter,...
                'Position',[160 360-dx*i wx wy].*WS,...
                'Style','edit');
            set(h2(i),'String','0');
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
                    
                    delete([h1(:); h2(:)])
                    spm_input(' ',1,'d')
                    break
                    
                end
            end
        end
    end
    delete(d)
    con=c(:);
    
otherwise,
    disp('Error in spm_dcm_contrasts: contrast must be for A, B or C');
    return
end

if nargin < 2
    close 
end


