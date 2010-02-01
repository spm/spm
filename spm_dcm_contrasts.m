function [con] = spm_dcm_contrasts(DCM_filename,D)
% Make contrast vector for a DCM model
% FORMAT [con] = spm_dcm_contrasts(DCM_filename,D)
%
% DCM_filename  DCM file name
% D             'A','B' or 'C' i.e. connectivity matrix of interest
%
% con           Column vector specifying contrast weights
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Will Penny
% $Id: spm_dcm_contrasts.m 3705 2010-02-01 20:51:28Z karl $
 
% Set-up
%--------------------------------------------------------------------------
Finter  = spm_figure('GetWin','Interactive');
header  = get(Finter,'Name');
set(Finter,'Name','Dynamic Causal Modelling')
WS      = spm('WinScale');
 
% load DCM if necessary
%--------------------------------------------------------------------------
try
    P   = DCM_filename;
    load(P);
catch
    DCM = DCM_filename;
end
 
Y.name  = DCM.Y.name;
U.name  = DCM.U.name;
Ep      = DCM.Ep;           % MAP estimates
n       = DCM.n;            % number of regions
m       = length(U.name);   % number of inputs
a       = zeros(n,n);
c       = zeros(n,m);
b       = zeros(n,n,m);
 
 
% prompt for contrast if necessary
%--------------------------------------------------------------------------
if (nargin < 2) | isempty(D)
    str = 'contrast for';
    D   = spm_input(str,1,'b',{'A','B','C'});
end
 
% set object sizes
%--------------------------------------------------------------------------
dx    = 35;
wx    = 30;
wy    = 20;
d     = uicontrol(Finter,'String','done','Position',[300 50 060 020].*WS);
 
% Define left edge dialogue
% itext_left = 080;
% inum_left  = 180;
%--------------------------------------------------------------------------
BackgroundColor = get(get(d,'Parent'),'Color');
itext_left  = 030;
inum_left   = 80;
name_length = 8; % number of characters in short names
for i=1:n
    if length(Y.name{i}) > name_length
        short_name(i).str = Y.name{i}(1:name_length);
    else
        short_name(i).str = Y.name{i};
    end
end
text_top = 336;
 
switch D
 
    % Get contrast weights
    %======================================================================
    case 'A' % intrinsic
 
        str   = sprintf('Enter contrast for A: ');
        spm_input(str,1,'d');
 
        % Print names and numbers of regions
        %------------------------------------------------------------------
        for i = 1:n
            str    = sprintf('%s   %i',short_name(i).str,i);
            h1(i)  = uicontrol(Finter,'String',str,...
                'Style','text',...
                'HorizontalAlignment','right',...
                'BackgroundColor',BackgroundColor,...
                'Position',[itext_left text_top-dx*i 080 020].*WS);
            h2(i)  = uicontrol(Finter,'String',sprintf('%i',i),...
                'Style','text',...
                'BackgroundColor',BackgroundColor,...
                'Position',[inum_left+dx*i text_top 020 020].*WS);
        end
 
        % Set contrast values to zero and display
        %------------------------------------------------------------------
        for i = 1:n
            for j = 1:n
                cc=ceil([inum_left+dx*j text_top+4-dx*i wx wy].*WS);
                h3(i,j) = uicontrol(Finter,...
                    'Position',cc,...
                    'BackgroundColor',BackgroundColor,...
                    'Style','edit');
                set(h3(i,j),'String','0');
 
            end
        end
        drawnow
 
        % wait for user to specify contrast weights and press 'done'
        %------------------------------------------------------------------
        while(1)
            pause(0.01)
            if strcmp(get(gco,'Type'),'uicontrol')
                if strcmp(get(gco,'String'),'done')
                    for i = 1:n
                        for j = 1:n
                            a(i,j) = str2num(get(h3(i,j),'string'));
                        end
                    end
                    delete([h1(:); h2(:); h3(:)])
                    spm_input(' ',1,'d')
                    break
                end
            end
        end
        con     = spm_unvec(spm_vec(Ep)*0,Ep);
        con.A   = a;
        con     = spm_vec(con);
 
    case 'B' % modulatory
        %------------------------------------------------------------------
        for k = 1:m,
            str   = sprintf(...
                'Enter contrast for B: effects of input  %-12s',...
                U.name{k});
            spm_input(str,1,'d')
 
            for i = 1:n
                str    = sprintf('%s   %i',short_name(i).str,i);
                h1(i)  = uicontrol(Finter,'String',str,...
                    'Style','text',...
                    'BackgroundColor',BackgroundColor,...
                    'Position',[itext_left text_top-dx*i 080 020].*WS);
                h2(i)  = uicontrol(Finter,'String',sprintf('%i',i),...
                    'Style','text',...
                    'BackgroundColor',BackgroundColor,...
                    'Position',[inum_left+dx*i text_top 020 020].*WS);
            end
            for i = 1:n
                for j = 1:n
                    cc=ceil([inum_left+dx*j text_top+4-dx*i wx wy].*WS);
                    h3(i,j) = uicontrol(Finter,...
                        'Position',cc,...
                        'BackgroundColor',BackgroundColor,...
                        'Style','edit');
                    set(h3(i,j),'String','0');
                end
            end
            drawnow
 
            % wait for 'done'
            %--------------------------------------------------------------
            set(gcf,'CurrentObject',h3(1))
            while(1)
                pause(0.01)
                if strcmp(get(gco,'Type'),'uicontrol')
                    if strcmp(get(gco,'String'),'done')
                        for i = 1:n
                            for j = 1:n
                                b(i,j,k) = str2num(get(h3(i,j),'string'));
                            end
                        end
                        delete([h1(:); h2(:); h3(:)])
                        spm_input(' ',1,'d')
                        break
 
                    end
                end
            end
 
        end
        con     = spm_unvec(spm_vec(Ep)*0,Ep);
        con.B   = b;
        con     = spm_vec(con);
 
    case 'C' % input
        %------------------------------------------------------------------
        for k = 1:m
            str   = sprintf(...
                'Enter contrast for C: Effects of input %-12s',...
                U.name{k});
            spm_input(str,1,'d');
 
            for i = 1:n
                str    = sprintf('%s   %i',short_name(i).str,i);
                h1(i)  = uicontrol(Finter,'String',str,...
                    'Style','text',...
                    'BackgroundColor',BackgroundColor,...
                    'Position',[itext_left text_top-dx*i 080 020].*WS);
                h2(i)  = uicontrol(Finter,...
                    'Position',[inum_left+dx text_top+4-dx*i wx wy].*WS,...
                    'BackgroundColor',BackgroundColor,...
                    'Style','edit');
                set(h2(i),'String','0');
            end
            drawnow
 
            % wait for 'done'
            %--------------------------------------------------------------
            set(gcf,'CurrentObject',h2(1))
            while(1)
                pause(0.01)
                if strcmp(get(gco,'Type'),'uicontrol')
                    if strcmp(get(gco,'String'),'done')
 
                        % get c
                        %--------------------------------------------------
                        for i = 1:n
                            c(i,k)   = str2num(get(h2(i),'string'));
                        end
 
                        delete([h1(:); h2(:)])
                        spm_input(' ',1,'d')
                        break
 
                    end
                end
            end
        end
        con     = spm_unvec(spm_vec(Ep)*0,Ep);
        con.C   = c;
        con     = spm_vec(con);
 
    otherwise
        disp('Error in spm_dcm_contrasts: contrast must be for A, B or C');
        close
        return
end
delete(d);
