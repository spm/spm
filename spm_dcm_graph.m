function spm_dcm_graph(xY,A)
% Region and anatomical graph display
% FORMAT spm_dcm_graph(xY,A)
% xY    - cell of region structures (see spm_regions)
% A     - connections of weighted directed graph
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_graph.m 4099 2010-10-22 19:47:37Z karl $
 
 
% display parameters
%--------------------------------------------------------------------------
col   = {'b','g','r','c','m','y','k','w'};
 
% get dimensions, locations and names
%--------------------------------------------------------------------------
m     = size(xY,2);
L     = [];
S     = [];
for i = 1:m
    L       = [L xY(i).xyz];
    name{i} = xY(i).name(1:min(end,3));
    S       = [S xY(i).spec];
end


% Render graph in anatomical space
%==========================================================================
subplot(2,1,1);cla
set(gca,'position',[0 .5 1 .5])
options.query = [];
options.hfig  = gcf;
options.ParentAxes = gca;
options.markersize = 32;
options.meshsurf   = fullfile(spm('Dir'),'canonical','cortex_5124.surf.gii');
h  = spm_eeg_displayECD(L,[],8,name,options);
for i = 1:m
    set(h.handles.ht(i),'FontWeight','bold')
end
set(h.handles.mesh,'FaceAlpha',1/16);
if nargin < 2, return, end
 
 
% Connections
%--------------------------------------------------------------------------
W     = max(abs(A),abs(A'));
W     = W - diag(diag(W));
W     = 3*W/max(W(:));
W     = W.*(W > 1/128);
for i = 1:length(A)
    for j = (i + 1):length(A)
        if W(i,j)
            
            % associate colour with the strongest influence
            %--------------------------------------------------------------
            if abs(A(i,j)) > abs(A(j,i)), c = j; else, c = i; end
            line(L(1,[i j]),L(2,[i j]),L(3,[i j]),'Color',col{c},...
                'LineStyle','-',...
                'LineWidth',W(i,j));
        end
    end
end


% Render graph in functional space
%==========================================================================
 
% Multidimensional scaling (with the Weighted Graph Laplacian)
%--------------------------------------------------------------------------
D      = diag(sum(W));
G      = D - W;
[U V]  = eig(pinv(G));
U      = U*sqrt(V);
[V i]  = sort(-diag(V));
U      = U(:,i(1:3))';
 
% Procrustean transform to align with anatomy (not currently used)
%--------------------------------------------------------------------------
U       = spm_detrend(U')';
L       = spm_detrend(L')';
% [R V S] = spm_svd(U*L');
% R       = R*S;
% U       = R'*U;
U       = diag(sign(diag(U*L')))*U;
U       = real(U*80/max(abs(U(:))));
 
 
subplot(2,1,2);cla
set(gca,'position',[0 0 1 .5])
options.ParentAxes = gca;
g      = spm_eeg_displayECD(U,[],16,name,options);
delete(g.handles.mesh)
delete(findobj(get(gcf,'Children'),'Type','uicontrol'))
for i = 1:m
    set(g.handles.ht(i),'FontWeight','bold')
end
 
% Connections
%--------------------------------------------------------------------------
for i = 1:length(A)
    for j = (i + 1):length(A)
        if W(i,j)
            
            % associate colour with the strongest influence
            %--------------------------------------------------------------
            if abs(A(i,j)) > abs(A(j,i)), c = j; else, c = i; end
            line(U(1,[i j]),U(2,[i j]),U(3,[i j]),'Color',col{c},...
                'LineStyle','-',...
                'LineWidth',W(i,j));
        end
    end
end

