function spm_MDP_plot(MDP)
% creates a movie of hierarchical expectations
% FORMAT spm_MDP_plot(MDP))
%
% DEM - {DEM} structures from visual search simulations
%     - (requires fields to specify the labels of states and outcomes)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_plot.m 6957 2016-12-01 10:37:00Z karl $

% Preliminaries
%--------------------------------------------------------------------------
clf; subplot(2,1,1),axis([2 20 0 8]);
title('Narrative construction','FontSize',16)
scale = 1 - 1/8;
set(gca,'XTick',[],'YTick',[]), box on

MDP.label.k       = 1;
MDP.label.factor  = {'story','where','decision'};
MDP.label.name{1} = {'flee wait feed wait',
    'wait wait wait feed'
    'wait flee wait feed'
    'flee wait feed flee'
    'wait wait wait flee'
    'wait flee wait flee'};

MDP.MDP.label.k          = 1;
MDP.MDP.label.factor     = {'what','where','flip','flip'};
MDP.MDP.label.name{1}    = {'flee','feed','wait'};
MDP.MDP.label.modality   = {'what','where'};
MDP.MDP.label.outcome{1} = {'null','bird','seed','cat'};
MDP.MDP.label.outcome{2} = {'1','2','3','4'};


% draw stuff
%--------------------------------------------------------------------------
try k1 = MDP.label.k; catch, k1 = 1; end
T1     = size(MDP.xn{k1},4);
F1     = fix(1024/(T1*numel(MDP.label.name{k1}{1})));
F1     = min(max(F1,8),32)
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
                h1(i,j) = text(16*i/nt,j/ns + 6,str(j),'Color',[1 1 1]*p(j,i),'FontSize',F1);
            end
        end
    end
    
    
    % next level down
    %----------------------------------------------------------------------
    if isfield(MDP,'mdp')
        try k2 = MDP.MDP.label.k; catch, k2 = 1; end
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
                        h2(i,j) = text(16*i/nt,j/ns + 4,str(j),'Color',[1 1 1]*p(j,i),'FontSize',F2);
                    end
                end
            end
            
            % next level down
            %--------------------------------------------------------------
            if isfield(MDP.mdp(t1),'kdem')
                
                
            else
                
                % observed outcomes
                %----------------------------------------------------------
                o      = MDP.mdp(t1).o;
                no     = size(o,1);
                T      = 16;
                F3     = fix(1024/(no*numel(MDP.MDP.label.outcome{1}{1})));
                F3     = min(max(F3,8),48);

                for t3 = 1:T
                    
                    % draw stuff
                    %------------------------------------------------------
                    for i = 1:no
                        str = MDP.MDP.label.outcome{i}{o(i,t2)};
                        try
                            set(h3(i,j),'String',str);
                        catch
                            h3(i,j) = text(16*i/no,2,str,'Color','r','FontSize',F3);
                        end
                    end
                    
                    % plot timelines
                    %------------------------------------------------------
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
                    
                    
                    % and save graphics
                    %------------------------------------------------------
                    try
                        M(end + 1) = getframe(gca);
                    catch
                        M = getframe(gca);
                    end
                    
                    
                end
            end
        end
    end
end

% save movie
%--------------------------------------------------------------------------
set(gca,'Userdata',{M,16})
set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')

