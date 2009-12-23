function [con_vec,con_mat] = spm_dcm_contrasts (DCM_filename,D)
% Make contrast vector for a DCM model 
% FORMAT [con_vec,con_mat] = spm_dcm_contrasts (DCM_filename,D)
%
% DCM_filename  DCM file name
% D             'A','B' or 'C' ie. contrast for which connectivity matrix
%
% con_vec       Column vector specifying contrast of parameters
% con_mat       The same contrast but in matrix format
%
%               The contrasts are also saved in the DCM structure
%
%               DCM.contrast().con_vec
%               DCM.contrast().con_mat
%               DCM.contrast().con_type ('A', 'B' or 'C')
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_dcm_contrasts.m 3654 2009-12-23 20:09:54Z karl $


Finter = spm_figure('GetWin','Interactive');
header = get(Finter,'Name');
set(Finter,'Name','Dynamic Causal Modelling')
WS     = spm('WinScale');

P = DCM_filename;
load(P);

Y.name   = DCM.Y.name;
U.name   = DCM.U.name;
n        = DCM.n;            % number of regions
m        = length(U.name);   % number of inputs
a        = zeros(n,n);
c        = zeros(n,m);
b        = zeros(n,n,m);

if (nargin < 2) | isempty(D)
    str     = 'contrast for';
    D       = spm_input(str,1,'b',{'A','B','C'});
end


dx    = 35;
wx    = 30;
wy    = 20;
d     = uicontrol(Finter,'String','done','Position',[300 50 060 020].*WS);

% Define left edge dialogue
% itext_left=080;
% inum_left=180;

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
text_top=336;

switch D,
    
case 'A' % intrinsic
    
    str   = sprintf(...
            'Enter contrast for A: ');
        spm_input(str,1,'d');
        
    % Print names and numbers of regions
    for i = 1:n
        str    = sprintf('%s   %i',short_name(i).str,i);
        h1(i)  = uicontrol(Finter,'String',str,...
            'Style','text',...
            'HorizontalAlignment','right',...
            'Position',[itext_left text_top-dx*i 080 020].*WS);
        h2(i)  = uicontrol(Finter,'String',sprintf('%i',i),...
            'Style','text',...
            'Position',[inum_left+dx*i text_top 020 020].*WS);
    end
    
    % Set contrast values to zero and display
    for i = 1:n
        for j = 1:n
            cc=ceil([inum_left+dx*j text_top+4-dx*i wx wy].*WS);
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
    con_mat = a;
    
case 'B' % modulatory
    %---------------------------------------------------
    
    for k = 1:m,
        str   = sprintf(...
            'Enter contrast for B: effects of input  %-12s',...
            U.name{k});
        spm_input(str,1,'d')
        
        for i = 1:n
            str    = sprintf('%s   %i',short_name(i).str,i);
            h1(i)  = uicontrol(Finter,'String',str,...
                'Style','text',...
                'Position',[itext_left text_top-dx*i 080 020].*WS);
            h2(i)  = uicontrol(Finter,'String',sprintf('%i',i),...
            'Style','text',...
            'Position',[inum_left+dx*i text_top 020 020].*WS);
        end
        for i = 1:n
            for j = 1:n
                cc=ceil([inum_left+dx*j text_top+4-dx*i wx wy].*WS);
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
    con_mat = b;
    
case 'C' % input
    %---------------------------------------------------
    
    for k = 1:m,
        str   = sprintf(...
            'Enter contrast for C: Effects of input %-12s',...
            U.name{k});
        spm_input(str,1,'d');
        
        for i = 1:n
            str    = sprintf('%s   %i',short_name(i).str,i);
            h1(i)  = uicontrol(Finter,'String',str,...
                'Style','text',...
                'Position',[itext_left text_top-dx*i 080 020].*WS);
            h2(i)  = uicontrol(Finter,...
                'Position',[inum_left+dx text_top+4-dx*i wx wy].*WS,...
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
    con_mat=c;
    
    
otherwise,
    disp('Error in spm_dcm_contrasts: contrast must be for A, B or C');
    close
    return
end

delete(d);

% disp('Contrast:');
% disp(con_mat);
con_vec=con_mat(:);


% Save contrast in DCM file
if isfield(DCM,'contrast')
    num_contrast=length(DCM.contrast)+1;
else
    num_contrast=1;
end
DCM.contrast(num_contrast).con_vec  = con_vec;
DCM.contrast(num_contrast).con_mat  = con_mat;
DCM.contrast(num_contrast).con_type = D;

if spm_matlab_version_chk('7') >= 0
    save(P(:),'-V6','DCM');
else
    save(P(:),'DCM');
end;

if nargin < 2
    close 
end


