function spm_MDP_factor_graph(MDP)
% Draws a factor graph corresponding to MDP
% FORMAT spm_MDP_factor_graph(MDP)
%
% MDP.A{G}(O,N1,...,NF) - likelihood of O outcomes given hidden states
% MDP.B{F}(NF,NF,MF)    - transitions among states under MF control states
%
% This routine draws a simplified (normal style) factor graph  based upon
% the size of likelihood and prior probability matrices (and labels). The
% resulting graph can either be interpreted in terms of a factor graph with
% factors corresponding to white boxes. Alternatively, it can be
% interpreted as a graphical model with coloured boxes corresponding to
% random variables. The magenta boxes denote outcomes (at intermediate
% levels of deep models, if specified). The cyan boxes denote hidden states
% and the puce boxes represent policies. If a hidden state is controllable
% (i.e., has more than one control dependent probability transition matrix)
% it is labelled in blue (and the hidden states are shown as a stack of
% boxes to indicate they are conditioned on several policies). Key message
% passing is illustrated with three arrows, corresponding to ascending
% likelihoods, forward and backward messages (1, 2 and 3,respectively).
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_factor_graph.m 7651 2019-08-03 12:35:15Z karl $


% deal with a sequence of trials
%==========================================================================
spm_figure('GetWin','Factor Graph'); clf


% options
%--------------------------------------------------------------------------
MDP = spm_MDP_check(MDP);

if isfield(MDP,'link')
    
    % deep model
    %----------------------------------------------------------------------
    [hg,af] = spm_MDP_factor_level(MDP,1/2);
    [ag,hf] = spm_MDP_factor_level(MDP.MDP,1/8);
    
    % link if levels
    %----------------------------------------------------------------------
    for g = 1:size(MDP.link,1)
        for f = 1:size(MDP.link,2)
            if ~isempty(MDP.link{g,f})
                pf = get(hg(f),'Position');
                pg = get(hf(g),'Position');
                annotation('line',[pg(1)+pg(3)/2 pf(1)+pf(3)/2],[pg(2)+pg(4) pf(2)]);
            end
        end
    end  
    
else
    
    % single level
    %----------------------------------------------------------------------
    spm_MDP_factor_level(MDP);
end

function [ag,af] = spm_MDP_factor_level(MDP,DY)
% FORMAT spm_MDP_factor_level(MDP,DY)
% draw factor graph for a single level
% MDP - structure
% DY  - vertical offset in normalised units (0,1)
%
% ag - handles to outcome boxes
% af - handles to factor boxes

% defaults
%--------------------------------------------------------------------------
if nargin < 2,  DY = 1/4; end

% numbers of transitions, policies and states
%--------------------------------------------------------------------------
Nf  = numel(MDP.B);                 % number of hidden state factors
Ng  = numel(MDP.A);                 % number of outcome factors
for f = 1:Nf
    Nu(f)    = size(MDP.B{f},3);    % number of hidden controls
    Ns(f)    = size(MDP.B{f},1);    % number of hidden states
end
for g = 1:Ng
    No(g)    = size(MDP.A{g},1);    % number of outcomes
end


% aspect ratio
%--------------------------------------------------------------------------
AR  = get(gcf,'Position');
AR  = AR(4)/AR(3);


%% draw factor graph
%--------------------------------------------------------------------------
dx    = 1/32;                       % width of boxes
dy    = 1/32/AR;                    % height of boxes
for t = fliplr(1:2)
    
    % outcome modalities
    %======================================================================
    for g = fliplr(1:Ng)
        
        % position
        %------------------------------------------------------------------
        x     = t*1/4 + (1 + 1/8)*g*dx - 1/8;
        y     = DY + (1 + 1/8/AR)*g*dy;
        
        % outcomes and names
        %------------------------------------------------------------------
        pos   = [x y dx dy];
        col   = [1 (1 - exp(-2*g/Ng)) 1];
        ag(g) = annotation('rectangle',pos,'FaceColor',col);
        if t == 2
            annotation('textbox',[x+1/32 y 1/8 dy],...
                'String',MDP.label.modality{g},...
                'FontWeight','bold',...
                'LineStyle','none',...
                'FontSize',10,'Color','m');
            str = MDP.label.outcome{g};
            if numel(str) > 4
                str = str(1:4);
                str{end + 1} = '...';
            end
            annotation('textbox',[x+1/32+1/16 y-1/64 1/8 dy],...
                'String',str,...
                'FontWeight','Normal',...
                'LineStyle','none',...
                'FontSize',8,'Color','m');
            
        end
        pos = [x y+1/16 dx dy];
        annotation('rectangle',pos,'FaceColor','w');
        annotation('textbox',pos,...
            'HorizontalAlignment','Center',...
            'VerticalAlignment','bottom',...
            'FontSize',12,...
            'FontWeight','bold',...
            'String','A');
        
        % edges
        %------------------------------------------------------------------
        annotation('line',[x+dx/2 x+dx/2],[y+dy y+1/16]);
        annotation('line',[x+dx/2 x+dx/2],[y+dy+1/16 y+1/16+1/32]);
        
    end
    
    for f = fliplr(1:Nf)
        
        % position
        %------------------------------------------------------------------
        x     = t*1/4 + (1 + 1/8)*f*dx  - 1/8;
        y     = DY + 1/8 + (1 + 1/8/AR)*f*dy;
        pos   = [x y dx dy];
        col   = [(1 - exp(-2*f/Nf)) 1 1];
        
        % factors and names
        %------------------------------------------------------------------
        for u = fliplr(1:Nu(f))
            annotation('rectangle',[x+(u - 1)*dx/8 y+(u - 1)*dy/8/AR dx dy],...
                'FaceColor',col);
        end
        
        % hidden states also for constraint
        %------------------------------------------------------------------
        annotation('rectangle',pos,'FaceColor',col);
        af(f) = annotation('textbox',pos,...
            'HorizontalAlignment','Center',...
            'VerticalAlignment','bottom',...
            'FontSize',12,...
            'FontWeight','bold',...
            'String','=');
        
        % names
        %------------------------------------------------------------------
        if t == 2
            h = annotation('textbox',[x+3/16 y+dy/2 1/8 dy],...
                'String',MDP.label.factor{f},...
                'FontWeight','bold',...
                'LineStyle','none',...
                'FontSize',10,'Color','c');
            if Nu(f) > 1
                set(h,'Color','b');
            end
            
            str = MDP.label.name{f};
            if numel(str) > 4
                str = str(1:4);
                str{end + 1} = '...';
            end
            h = annotation('textbox',[x+3/16+1/16 y+dy/2-1/64 1/8 dy],...
                'String',str,...
                'FontWeight','Normal',...
                'LineStyle','none',...
                'FontSize',8,'Color','c');
            if Nu(f) > 1
                set(h,'Color','b');
            end
            
        end
        if t == 1
            
            if f ==1
                annotation('textarrow',[x x+1/32]+1/16,[y y] + 3*dy/4,...
                    'String','2',...
                    'FontWeight','bold',...
                    'HeadStyle','plain',...
                    'LineWidth',2,...
                    'HeadWidth',8,...
                    'HeadLength',8,...
                    'FontSize',10,'Color','r');
                annotation('textarrow',[x+1/32 x]+3/16,[y y] + 3*dy/4,...
                    'String','3',...
                    'FontWeight','bold',...
                    'HeadStyle','plain',...
                    'LineWidth',2,...
                    'HeadWidth',8,...
                    'HeadLength',8,...
                    'FontSize',10,'Color','r');
                 annotation('textarrow',[x x]+dx/4,[y y+1/48/AR] - 1/32/AR,...
                    'String','1',...
                    'FontWeight','bold',...
                    'HeadStyle','plain',...
                    'LineWidth',2,...
                    'HeadWidth',8,...
                    'HeadLength',8,...
                    'FontSize',10,'Color','r');
                
            end
        end
        
        % model parameters of factor
        %------------------------------------------------------------------
        pos  = [x+1/8 y dx dy];
        annotation('rectangle',pos,'FaceColor','w');
        annotation('textbox',pos,...
            'HorizontalAlignment','Center',...
            'VerticalAlignment','bottom',...
            'FontSize',12,...
            'FontWeight','bold',...
            'String','B');
        
        % edges
        %------------------------------------------------------------------
        try
            annotation('line',[x+dx x+1/8],[y+dy/2 y+dy/2]);
            annotation('line',[x+dx/2 x+dx/2],[y y-1/32]);
        end
        try
            annotation('line',[x+1/8+dx x+1/4],[y+dy/2 y+dy/2]);
        end
        
        % policies and free energies
        %------------------------------------------------------------------
        if t < 2 && Nu(f) > 1
            pos  = [x+1/8 y+1/16 dx dy];
            annotation('rectangle',pos,'FaceColor','w');
            annotation('textbox',pos,...
                'HorizontalAlignment','Center',...
                'VerticalAlignment','bottom',...
                'FontSize',12,...
                'FontWeight','bold',...
                'String','G');
            pos  = [x+1/8 y+1/8 dx dy];
            col  = [(1 - exp(-2*f/Nf)) (1 - exp(-2*f/Nf)) 1];
            annotation('rectangle',pos,'FaceColor',col);
            annotation('textbox',pos,...
                'HorizontalAlignment','Center',...
                'VerticalAlignment','bottom',...
                'FontSize',12,...
                'FontWeight','bold',...
                'String','=');
            
            % edges
            %--------------------------------------------------------------
            annotation('line',[x+1/8+dx/2 x+1/8+dx/2],[y+1/16+dy y+1/8]);
            annotation('line',[x+dx/2 x+1/8],[y+1/16+dy/2 y+1/16+dy/2]);
            annotation('line',[x+1/8+dx x+1/4+dx/2],[y+1/16+dy/2 y+1/16+dy/2]);
            annotation('line',[x+dx/2 x+dx/2],[y+dy y+1/16+dy/2]);
            annotation('line',[x+1/4+dx/2 x+1/4+dx/2],[y+dy y+1/16+dy/2]);
            
        end
        
        
        % edges: from hidden states to outcomes
        %------------------------------------------------------------------
        x1     = t*1/4 + (1 + 1/8)*dx - 1/8;
        y1     = DY + (1 + 1/8/AR)*dy;
        xg     = t*1/4 + (1 + 1/8)*max(Nf,Ng)*dx - 1/8;
        yg     = DY + (1 + 1/8/AR)*max(Nf,Ng)*dy;
        
        annotation('line',[x1,xg] + dx/2,[y1,yg] + 1/16 + 1/32);
        
    end
end
