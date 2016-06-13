function con = spm_dcm_connectivity_ui(DCM,D,title_text,defaults,enabled)
% GUI for manually specifying connection values in a DCM
% FORMAT [con] = spm_dcm_connectivity_ui(DCM,D,title_text,defaults,enabled)
%
% DCM        - DCM structure
% D          - 'A','B' or 'C' i.e. connectivity matrix of interest
% title_text - Text to display above the matrix, e.g. 'Enter contrast for '
% defaults   - (optional) structure of default values containing
%              defaults.A, defaults.B and defaults.C
% enabled    - (optional) structure of inputs to enable with binary 
%              matrices enabled.A, enabled.B and enabled.C
% 
% Returns:
% con        - structure with con.A, con.B and con.C of user-entered values
%__________________________________________________________________________
% Copyright (C) 2014-2016 Wellcome Trust Centre for Neuroimaging    

% Will Penny & Peter Zeidman
% $Id: spm_dcm_connectivity_ui.m 6808 2016-06-13 16:48:30Z guillaume $


% Set-up data
%--------------------------------------------------------------------------
n       = DCM.n;              % number of regions
m       = length(DCM.U.name); % number of inputs

region_names = DCM.Y.name;
input_names  = DCM.U.name;

try, A = defaults.A; catch, A = zeros(n,n); end
try, B = defaults.B; catch, B = zeros(n,n,m); end
try, C = defaults.C; catch, C = zeros(n,m); end

if nargin < 5
    enabled = struct('A',ones(size(A)),'B',ones(size(B)),'C',ones(size(C)));
end

% Set-up GUI
%--------------------------------------------------------------------------
Finter  = spm_figure('GetWin','Interactive');
header  = get(Finter,'Name');
set(Finter,'Name','Dynamic Causal Modelling')
WS      = spm('WinScale');
 
% prompt for contrast if necessary
%--------------------------------------------------------------------------
if (nargin < 2) || isempty(D)
    str = 'contrast for';
    D   = spm_input(str,1,'b',{'A','B','C'});
end
 
% set object sizes
%--------------------------------------------------------------------------
dx    = 35;
wx    = 30;
wy    = 20;
uic   = uicontrol(Finter,'String','done','Position',[300 50 060 020].*WS,...
        'Callback', 'uiresume(gcbf)');
 
% Define left edge dialogue
% itext_left = 080;
% inum_left  = 180;
%--------------------------------------------------------------------------
BackgroundColor = get(get(uic,'Parent'),'Color');
itext_left  = 030;
inum_left   = 80;
name_length = 8; % number of characters in short names
for i=1:n
    if length(region_names{i}) > name_length
        short_name(i).str = region_names{i}(1:name_length);
    else
        short_name(i).str = region_names{i};
    end
end
text_top = 336;

switch D
 
    % Get contrast weights
    %======================================================================
    case 'A' % intrinsic
 
        str   = sprintf('%s A: ',title_text);
        spm_input(str,1,'d');
        
        % Warn if no matrix elements enabled
        any_enabled = any(enabled.A(:));
        if any_enabled
            h0 = [];
        else
            h0 = add_none_enabled_warning();
        end        
 
        % Print names and numbers of regions
        %------------------------------------------------------------------
        for i = 1:n
            str    = sprintf('%s   %i',short_name(i).str,i);
            
            h1(i) = add_text([itext_left text_top-dx*i 080 020].*WS,...
                        str, any_enabled, 'right');
                    
            h2(i) = add_text([inum_left+dx*i text_top 020 020].*WS,...
                        sprintf('%i',i), any_enabled);
        end
 
        % Set contrast values and display
        %------------------------------------------------------------------
        for i = 1:n
            for j = 1:n
                cc=ceil([inum_left+dx*j text_top+4-dx*i wx wy].*WS);                          
                h3(i,j) = add_matrix_input(cc, num2str(A(i,j)), enabled.A(i,j)==1);
            end
        end
        drawnow
 
        % wait for 'done'
        %------------------------------------------------------------------
        uiwait(Finter);        
        for i = 1:n
            for j = 1:n
                A(i,j) = get_matrix_input(h3(i,j));
            end
        end
        
        % clean up
        %------------------------------------------------------------------        
        delete([h0; h1(:); h2(:); h3(:)]);
        spm_input(' ',1,'d');
         
    case 'B' % modulatory
        %------------------------------------------------------------------
        for k = 1:m,
            str   = sprintf(...
                '%s B: effects of input  %-12s',...
                title_text, input_names{k});
            spm_input(str,1,'d')
 
            % Warn if no matrix elements enabled
            any_enabled = any(any(enabled.B(:,:,k)));
            if any_enabled
                h0 = [];
            else
                h0 = add_none_enabled_warning();
            end
            
            % Print names and numbers of regions
            %--------------------------------------------------------------
            for i = 1:n
                str    = sprintf('%s   %i',short_name(i).str,i);
       
                h1(i) = add_text([itext_left text_top-dx*i 080 020].*WS,...
                        str, any_enabled);                
                    
                h2(i) = add_text([inum_left+dx*i text_top 020 020].*WS,...
                        sprintf('%i',i), any_enabled);                       
            end
            
            % Set contrast values and display
            %--------------------------------------------------------------
            for i = 1:n
                for j = 1:n
                    cc=ceil([inum_left+dx*j text_top+4-dx*i wx wy].*WS);
                    h3(i,j) = add_matrix_input(cc, num2str(B(i,j,k)), enabled.B(i,j,k)==1);                  
                end
            end
            drawnow
 
            % wait for 'done'
            %--------------------------------------------------------------
            set(gcf,'CurrentObject',h3(1));
            uiwait(Finter);   
            for i = 1:n
                for j = 1:n
                    B(i,j,k) = get_matrix_input(h3(i,j));
                end
            end
            
            % clean up
            %--------------------------------------------------------------
            delete([h0; h1(:); h2(:); h3(:)]);
            spm_input(' ',1,'d');
 
        end
 
    case 'C' % input
        %------------------------------------------------------------------
        for k = 1:m
            str   = sprintf(...
                '%s C: Effects of input %-12s',...
                title_text, input_names{k});
            spm_input(str,1,'d');
            
            % Warn if no matrix elements enabled
            any_enabled = any(enabled.C(:,k));
            if any_enabled
                h0 = [];
            else
                h0 = add_none_enabled_warning();
            end                                   
            
            % Print names and numbers of regions
            %--------------------------------------------------------------
            for i = 1:n
                str    = sprintf('%s   %i',short_name(i).str,i);

                h1(i) = add_text([itext_left text_top-dx*i 080 020].*WS,...
                        str, any_enabled);                      

                cc = [inum_left+dx text_top+4-dx*i wx wy].*WS;
                h2(i) = add_matrix_input(cc, num2str(C(i,k)), enabled.C(i,k)==1);                  

            end
            drawnow
 
            % wait for 'done'
            %--------------------------------------------------------------
            set(gcf,'CurrentObject',h2(1))
            uiwait(Finter);   

            for i = 1:n
                C(i,k)   = get_matrix_input(h2(i));
            end

            % clean up
            %--------------------------------------------------------------            
            delete([h0; h1(:); h2(:)]);
            spm_input(' ',1,'d');
        end
    otherwise
        disp('Error in spm_dcm_connectivity_ui: matrix must be A, B or C');
        close
        return
end

% Bundle outputs
%--------------------------------------------------------------------------
con = struct('A',A,'B',B,'C',C);

delete(uic);

%--------------------------------------------------------------------------
function h = add_matrix_input(position, default_txt, enabled)
    % Adds an editor for a single element of a matrix
    if enabled == 1
        isvisible = 'On';
    else
        isvisible = 'Off';
    end            
    
    h = uicontrol(Finter,...
        'Position',position,...
        'BackgroundColor',BackgroundColor,...
        'Style','edit',...
        'Visible',isvisible,...
        'String',default_txt);
end
%--------------------------------------------------------------------------
function num = get_matrix_input(h)
    % Reads a number from a single element of a matrix
    entry=get(h,'string');
    entry=strtrim(entry);
    
    if isempty(entry)
        entry = '0';
    end
    
    try
        num = str2double(entry);
    catch
        error('Could not read matrix input');
    end
end
%--------------------------------------------------------------------------
function h = add_text(position, default_txt, visible, horizontalAlignment)
    % Adds a text label to the GUI
    if nargin < 3
        visible = 1;
    end
    if nargin < 4
        horizontalAlignment = 'left';
    end
    
    if visible
        visible = 'On';
    else
        visible = 'Off';
    end
    
    h  = uicontrol(Finter,'String',default_txt,...
        'Style','text',...
        'HorizontalAlignment',horizontalAlignment,...
        'BackgroundColor',BackgroundColor,...
        'Position',position,...
        'Visible',visible);
end
%--------------------------------------------------------------------------
function h = add_none_enabled_warning()
    % Adds a warning for when there are no enabled matrix entries
    position = [itext_left text_top-dx*1 250 020].*WS;
    nonestring = 'Nothing to enter - press done to continue';
    h = add_text(position, nonestring);
end

end
