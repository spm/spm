function spm_MDP_plot(MDP)
% creates a movie of hierarchical expectations and outcomes
% FORMAT spm_MDP_plot(MDP))
%
% MDP - nested MDP (and DEM) structures
%     - (requires fields to specify the labels of states and outcomes)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_plot.m 6977 2016-12-24 17:48:44Z karl $

% Preliminaries
%--------------------------------------------------------------------------
cla,axis([1 20 0 8]); axis image, hold on, box on
set(gca,'Fontsize',16,'XColor','r','TickLength',[0 0],'XTick',[],'YTick',[])
title('Narrative construction','FontSize',16)

% switch off basis functions for visual stimuli
%--------------------------------------------------------------------------
try
global STIM
B = STIM.B;
STIM.B = 1;
end

% hidden factors to display (one per level)
%--------------------------------------------------------------------------
scale  = 1 - 1/16;                             % background level for text
try k1 = MDP.label.k;     catch, k1 = 1; end
try k2 = MDP.MDP.label.k; catch, k2 = 1; end

% cycle over highest level
%--------------------------------------------------------------------------
T1     = size(MDP.xn{k1},4);
F1     = fix(1024/(T1*numel(MDP.label.name{k1}{1})));
F1     = min(max(F1,8),32);
for t1 = 1:T1
    
    % draw stuff
    %----------------------------------------------------------------------
    p     = 1 - squeeze(MDP.xn{k1}(end,:,:,t1));
    p     = p*scale;
    str   = MDP.label.name{k1};
    ns    = size(p,1);
    nt    = size(p,2);
    for i = 1:nt
        for j = 1:ns
            try
                set(h1(i,j),'Color',[1 1 1]*p(j,i));
            catch
                h1(i,j) = text(16*i/nt,j/ns + 6,str(j),...
                    'HorizontalAlignment','center',...
                    'Color',[1 1 1]*p(j,i),'FontSize',F1);
            end
        end
    end
    
    
    % next level down
    %----------------------------------------------------------------------
    if isfield(MDP,'mdp')
        T2     = size(MDP.mdp(t1).xn{k2},4);
        F2     = fix(1024/(T2*numel(MDP.MDP.label.name{k1}{1})));
        F2     = min(max(F2,8),32);
        for t2 = 1:T2
            
            % draw stuff
            %--------------------------------------------------------------
            p     = 1 - squeeze(MDP.mdp(t1).xn{k2}(end,:,:,t2));
            p     = p*scale;
            str   = MDP.MDP.label.name{k2};
            ns    = size(p,1);
            nt    = size(p,2);
            for i = 1:nt
                for j = 1:ns
                    try
                        set(h2(i,j),'Color',[1 1 1]*p(j,i));
                    catch
                        h2(i,j) = text(16*i/nt,j/ns + 4,str(j),...
                            'HorizontalAlignment','center',...
                            'Color',[1 1 1]*p(j,i),'FontSize',F2);
                    end
                end
            end
            
            % check for extra time points
            %--------------------------------------------------------------
            try
                delete(h2((nt + 1):end,:))
            end
            
            % next level down
            %--------------------------------------------------------------
            if isfield(MDP.mdp(t1),'dem')
                
                % accumulated expectations from Bayesian filtering
                %----------------------------------------------------------
                o  = MDP.mdp(t1).dem(t2).X;
                no = numel(o);
                T  = size(o{1},2);
                
            else
                
                % observed outcomes
                %----------------------------------------------------------
                o  = MDP.mdp(t1).o;
                no = size(o,1);
                T  = 16;
            end
            
            % outcome modalities
            %--------------------------------------------------------------
            set(gca,'XTickLabel',MDP.MDP.label.modality)
            set(gca,'XTick',16*(1:no)/no)
            
            F3     = fix(1024/(no*numel(MDP.MDP.label.outcome{1}{1})));
            F3     = min(max(F3,8),32);
            for t3 = 1:T
                
                % draw stuff
                %----------------------------------------------------------
                for i = 1:no
                    
                    if iscell(o)
                        
                        % probabilistic outcomes (expectations)
                        %--------------------------------------------------
                        x     = 16*i/no;
                        p     = spm_softmax(log(o{i}(:,t3))/4);
                        str   = MDP.MDP.label.outcome{i};
                        [p,j] = sort(p,'descend');
                        col   = [1 1 1] - [0 1 1]*(p(1) - 1/numel(p));
                        try
                            set(h3(i),'String',str(j(1)),...
                                'Color',col);
                        catch
                            h3(i) = text(x,3,str(j(1)),...
                                'HorizontalAlignment','center',...
                                'Color',col,'FontSize',F3);
                        end
                        
                        % observations in image format
                        %--------------------------------------------------
                        x     = [-1 1] + x;
                        y     = [-1 1] + 1;
                        try
                            px = MDP.mdp(t1).dem(t2).pU.x{1}(:,t3);
                            pv = MDP.mdp(t1).dem(t2).pU.v{2}(:,t3);
                            px = spm_unvec(px,MDP.MDP.DEM.G(1).x);
                            pv = spm_unvec(pv,MDP.MDP.DEM.G(2).v);
                            Y  = MDP.MDP.DEM.label{i}(px,pv);
                            try
                                delete(h4(i))
                            end
                            h4(i) = image(x,y,flipud(Y)*64);
                            
                        end
                        
                    else
                        
                        % deterministic outcomes
                        %--------------------------------------------------
                        str = MDP.MDP.label.outcome{i}{o(i,t2)};
                        try
                            set(h3(i),'String',str);
                        catch
                            h3(i) = text(16*i/no,2,str,'Color','r',...
                                'HorizontalAlignment','center','FontSize',F3);
                        end
                    end
                end
                
                % plot timelines
                %----------------------------------------------------------
                x  = [1 1]*(t3 + (t2 - 1)*T + (t1 - 1)*T*MDP.MDP.T);
                x  = 1 + 16*x/(MDP.MDP.T*MDP.T*T);
                y  = [7.2 8];
                try
                    set(l1,'XData',x);
                catch
                    l1 = line(x,y,'Color','r','LineWidth',8);
                end
                
                x  = [1 1]*(t3 + (t2 - 1)*T);
                x  = 1 + 16*x/(T*T2);
                y  = [5.2 6];
                try
                    set(l2,'XData',x);
                catch
                    l2 = line(x,y,'Color','r','LineWidth',8);
                end
                
                x  = [1 1]*t3;
                x  = 1 + 16*x/T;
                y  = [3.2 4];
                try
                    set(l3,'XData',x);
                catch
                    l3 = line(x,y,'Color','r','LineWidth',8);
                end
                
                % and save graphics
                %----------------------------------------------------------
                try
                    M(end + 1) = getframe(gca);
                catch
                    M = getframe(gca);
                end
                
            end
        end
        
    else
        
        % observed outcomes
        %------------------------------------------------------------------
        o      = MDP.o;
        no     = size(o,1);
        T      = 16;
        F3     = fix(1024/(no*numel(MDP.label.outcome{1}{1})));
        F3     = min(max(F3,8),48);
        
        % outcome modalities
        %------------------------------------------------------------------
        set(gca,'XTickLabel',MDP.label.modality)
        set(gca,'XTick',16*(1:no)/no)
        
        for t3 = 1:T
            
            % draw stuff
            %--------------------------------------------------------------
            for i = 1:no
                str = MDP.label.outcome{i}{o(i,t1)};
                try
                    set(h3(i,j),'String',str);
                catch
                    h3(i,j) = text(16*i/no,2,str,'Color','r',...
                        'HorizontalAlignment','center','FontSize',F3);
                end
            end
            
            % plot timelines
            %--------------------------------------------------------------
            x  = [1 1]*(t3 + (t1 - 1)*T);
            x  = 1 + 16*x/(T*T1);
            y  = [7.2 8];
            try
                set(l2,'XData',x);
            catch
                l2 = line(x,y,'Color','r','LineWidth',8);
            end
            
            x  = [1 1]*t3;
            x  = 1 + 16*x/T;
            y  = [5.2 6];
            try
                set(l3,'XData',x);
            catch
                l3 = line(x,y,'Color','r','LineWidth',8);
            end
            
            % and save graphics
            %--------------------------------------------------------------
            try
                M(end + 1) = getframe(gca);
            catch
                M = getframe(gca);
            end
            
        end
    end
end


% save movie
%--------------------------------------------------------------------------
try, STIM.B = B; end
set(gca,'Userdata',{M,16})
set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')

