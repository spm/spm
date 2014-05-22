function spm_dcm_graph(xY,A)
% Region and anatomical graph display
% FORMAT spm_dcm_graph(xY,[A])
% xY    - cell of region structures (see spm_regions) (fMRI)
%         or ECD locations xY.Lpos and xY.Sname (EEG)
% A     - connections of weighted directed graph
%__________________________________________________________________________
% Copyright (C) 2010-2014 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_graph.m 6012 2014-05-22 18:21:41Z guillaume $


% get dimensions, locations and names
%--------------------------------------------------------------------------
try
    % fMRI
    %----------------------------------------------------------------------
    m     = size(xY,2);
    L     = [];
    for i = 1:m
        L       = [L xY(i).xyz];
        name{i} = xY(i).name(1:min(end,3));
    end
    
catch
    
    % EEG
    %----------------------------------------------------------------------
    L    = xY.Lpos;
    name = xY.Sname;
    
end

% display parameters
%--------------------------------------------------------------------------
col   = {'b','g','r','c','m','y','k','w'};
m     = size(L,2);


%-Render graph in anatomical space
%==========================================================================
ax = subplot(2,1,1); cla(ax);
set(ax,'position',[0 .5 1 .5])
options.query = [];
options.hfig  = ancestor(ax,'figure');
options.ParentAxes = ax;
options.markersize = 32;
options.meshsurf = fullfile(spm('Dir'),'canonical','iskull_2562.surf.gii');
spm_eeg_displayECD(L(:,1),[],0,[],options);
options.meshsurf = fullfile(spm('Dir'),'canonical','cortex_8196.surf.gii');
h = spm_eeg_displayECD(L,[],8,name,options);
set(h.handles.ht,'FontWeight','bold')
set(h.handles.mesh,'FaceAlpha',1/16);


% return if no connectivity
%--------------------------------------------------------------------------
if nargin < 2, return, end

% check for cell array connections (EEG)
%--------------------------------------------------------------------------
if iscell(A)
    W = 0;
    C = 0;
    for i = 1:length(A)
        C = C + abs(exp(A{i}));
        W = W + max(C,C');
    end
    A = C;
    
elseif isnumeric(A)
    W = max(abs(A),abs(A'));
    
else
    W = [];
end

% Connections - if weights (W) are defined
%--------------------------------------------------------------------------
if numel(W)
    
    W     = W - diag(diag(W));
    W     = 3*W/max(W(:));
    W     = W.*(W > 1/128);
    for i = 1:length(A)
        for j = (i + 1):length(A)
            if W(i,j)
                
                % associate colour with the strongest influence
                %----------------------------------------------------------
                if abs(A(i,j)) > abs(A(j,i)), c = j; else c = i; end
                k   = rem(c - 1,length(col)) + 1;
                line(L(1,[i j]),L(2,[i j]),L(3,[i j]),'Color',col{k},...
                    'LineStyle','-',...
                    'LineWidth',W(i,j));
            end
        end
    end
end


%-Render graph in functional space (with the locations U)
%==========================================================================

if isstruct(A)
    
    P.A    = zeros(m,m);
    P      = spm_dcm_fmri_graph_gen([],A,P);
    W      = P.A;
    
    i      = 1:min(size(A.x,1),3);
    U      = zeros(3,m);
    U(i,:) = A.x(i,:);
    A      = W;
    W      = sign(W);
    
else
      
    % Multidimensional scaling (with the Weighted Graph Laplacian)
    %----------------------------------------------------------------------
    D      = diag(sum(W));
    G      = D - W;
    [U,V]  = eig(full(spm_pinv(G)));
    U      = U*sqrt(V);
    [V,i]  = sort(-diag(V));
    U      = U(:,i(1:3))';
    
end

% Procrustean transform
%----------------------------------------------------------------------
U      = spm_detrend(U')';
U      = real(U*40/max(abs(U(:))));


ax = subplot(2,1,2); cla(ax);
set(ax,'position',[0 0 1 .5])
options.ParentAxes = ax;
if m > 8; i = 8; else i = 16; end
g     = spm_eeg_displayECD(U,[],i,name,options);
delete(g.handles.mesh)
delete(findobj(get(gcf,'Children'),'Type','uicontrol'))
for i = 1:m
    set(g.handles.ht(i),'FontWeight','bold')
end

% Connections
%--------------------------------------------------------------------------
for i = 1:m
    for j = (i + 1):m
        
        % associate colour with the strongest influence
        %------------------------------------------------------------------
        if abs(A(i,j)) > abs(A(j,i)), c = j; else c = i; end
        k   = rem(c - 1,length(col)) + 1;
        
        if W(i,j) > 0
            line(U(1,[i j]),U(2,[i j]),U(3,[i j]),'Color',col{k},...
                'LineStyle','-',...
                'LineWidth', W(i,j));
        elseif W(i,j) < 0
            line(U(1,[i j]),U(2,[i j]),U(3,[i j]),'Color',col{k},...
                'LineStyle','-.',...
                'LineWidth',-W(i,j));  
        end
    end
end
