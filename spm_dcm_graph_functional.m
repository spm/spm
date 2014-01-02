function spm_dcm_graph_functional(A,V)
% Functional graph display
% FORMAT spm_dcm_graph_functional(A,V)
% FORMAT spm_dcm_graph_functional(V) - metric MDS
% A     - (m x m) weighted adjacency matrix
% V     - (n x m) locations in (nD) Multidimensional Scaling (MDS) Space 
%
% If V is not specified the Weighted Graph Laplacian of A is used with
% metric MDS to define the functional space.
%__________________________________________________________________________
% Copyright (C) 2013-2014 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_graph_functional.m 5823 2014-01-02 14:01:10Z guillaume $


% colours and number of nodes
%--------------------------------------------------------------------------
col   = {'b','g','r','c','m','y','k','w'};
m     = size(A,2);

% metric MDS
%--------------------------------------------------------------------------
if nargin == 1
    
    % unit sphere
    %----------------------------------------------------------------------
    set(gcf,'renderer','OpenGL');
    [x,y,z] = sphere(64); cla
    surf(x,y,z,'edgecolor','none','facealpha',1/8), hold on
    
    i      = 1:min(size(A,1),3);
    U      = zeros(3,m);
    U(i,:) = A(i,:);
    
    U      = U/sqrt(mean(sum(U.^2)));
    u      = spm_en(U);
    
    % nodes
    %----------------------------------------------------------------------
    for i = 1:m
        k   = rem(i - 1,length(col)) + 1;
        line(U(1,i),U(2,i),U(3,i),'Color',col{k},...
            'Marker','.',...
            'MarkerSize',32);
        line(u(1,i),u(2,i),u(3,i),'Color',col{k},...
            'Marker','.',...
            'MarkerSize',64);
        line([u(1,i) 0],[u(2,i) 0],[u(3,i) 0],'Color',col{k},...
            'LineStyle','-.');  
    end
    
    light
    lighting gouraud
    material shiny
    axis image
    return
end


% get wieghted (semi-defintie positive) adjacency matrix W from A
%==========================================================================

% check for cell array connecions (EEG)
%--------------------------------------------------------------------------
if iscell(A)
    W = 0;
    C = 0;
    for i = 1:length(A)
        C = C + abs(exp(A{i}));
        W = W + max(C,C');
    end
    A = C;
else
    try
        W = max(abs(A),abs(A'));
    end
end


% Connections - if weights (W) are defined
%--------------------------------------------------------------------------
if nargin < 2
    
    % Multidimensional scaling (with the Weighted Graph Laplacian)
    %----------------------------------------------------------------------
    W      = W - diag(diag(W));
    D      = diag(sum(W));
    G      = D - W;
    [U,V]  = eig(full(spm_pinv(G)));
    U      = U*sqrt(V);
    [V,i]  = sort(-diag(V));
    U      = U(:,i(1:3))';
    
else
    
    i      = 1:min(size(V,1),3);
    U      = zeros(3,m);
    U(i,:) = V(i,:);
    
end


%-Render graph in functional space (with the locations U)
%==========================================================================

% Procrustean transform
%--------------------------------------------------------------------------
U      = spm_detrend(U')';
U      = U/max(abs(U(:)));
A      = A/max(abs(A(:)));


% Nodes and connections
%--------------------------------------------------------------------------
cla;
for i = 1:m
    
    % nodes
    %----------------------------------------------------------------------
    k   = rem(i - 1,length(col)) + 1;
    line(U(1,i),U(2,i),U(3,i),'Color',col{k},...
        'LineStyle','.',...
        'MarkerSize',64);
    
    % edges
    %----------------------------------------------------------------------
    for j = (i + 1):m
        if A(i,j) > 1/16
            line(U(1,[i j]),U(2,[i j]),U(3,[i j]),'Color','k',...
                'LineStyle','-');
        elseif  A(i,j) < -1/16
            line(U(1,[i j]),U(2,[i j]),U(3,[i j]),'Color','k',...
                'LineStyle','-.');
        end
    end
end

axis image off
hold off
