function [DCM] = spm_dcm_phase_results(DCM,Action)
% Results for Dynamic Causal Modeling (DCM) for phase coupling
% FORMAT spm_dcm_phase_results(DCM,Action);
% Action:
%     'Sin(Data) - Region j'
%     'Coupling (As)'
%     'Coupling (Bs)'
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging
 
% Will Penny
% $Id: spm_dcm_phase_results.m 3472 2009-10-16 17:10:26Z will $


% get figure handle
%--------------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');
colormap(gray)
figure(Fgraph)
clf

xY     = DCM.xY;
nt     = length(xY.y);           % Nr of trial types
nr     = size(DCM.As,1);          % Nr of sources
nu     = length(DCM.B);          % Nr of experimental effects
ns     = size(xY.y{1},1);        % Nr of time bins
pst    = xY.pst;                 % peri-stmulus time

    
% switch
%--------------------------------------------------------------------------
switch(lower(Action))
    
    case{lower('Coupling (As)')}
    
   
    dx=min(0.1,1/nr);
    hold on
    
    cx=0;
    cy=0;
    text(cx,cy,'Posterior mean');
    % reconstitute time-frequency coupling
    %----------------------------------------------------------------------
    for i = 1:nr
        for j = 1:nr
            if (i==1)
                text(i*dx+cx,(j+1)*dx+cy,DCM.Sname{j});
            end
            if (j==1)
                text((i+1)*dx+cx,j*dx+cy,DCM.Sname{i});
            end
            text((i+1)*dx+cx,(j+1)*dx+cy,sprintf('%1.3f',abs(DCM.Ep.As(j,i))));    
        end
    end
    axis ij
    axis off
    title('Endogenous coupling (As)')
     
    cx=0;
    cy=0.5;
    
    text(cx,cy,'Posterior probability |As| > 0');
    for i = 1:nr
        for j = 1:nr
            if (i==1)
                text(i*dx+cx,(j+1)*dx+cy,DCM.Sname{j});
            end
            if (j==1)
                text((i+1)*dx+cx,j*dx+cy,DCM.Sname{i});
            end
            text((i+1)*dx+cx,(j+1)*dx+cy,sprintf('%1.3f',DCM.Pp.As(j,i)));    
        end
    end
    axis ij
    axis off
    
case{lower('Coupling (Bs)')}
    
    dx=min(0.1,1/nr);
    hold on
    
    cx=0;
    cy=0;
    text(cx,cy,'Posterior mean');
    % reconstitute time-frequency coupling
    %----------------------------------------------------------------------
    for i = 1:nr
        for j = 1:nr
            if (i==1)
                text(i*dx+cx,(j+1)*dx+cy,DCM.Sname{j});
            end
            if (j==1)
                text((i+1)*dx+cx,j*dx+cy,DCM.Sname{i});
            end
            text((i+1)*dx+cx,(j+1)*dx+cy,sprintf('%1.3f',DCM.Ep.Bs{1}(j,i)));    
        end
    end
    axis ij
    axis off
    title({'changes in coupling (Bs)';DCM.xU.name{1}})
    
    cx=0;
    cy=0.5;
    
    text(cx,cy,'Posterior probability |Bs| > 0');
    for i = 1:nr
        for j = 1:nr
            if (i==1)
                text(i*dx+cx,(j+1)*dx+cy,DCM.Sname{j});
            end
            if (j==1)
                text((i+1)*dx+cx,j*dx+cy,DCM.Sname{i});
            end
            text((i+1)*dx+cx,(j+1)*dx+cy,sprintf('%1.3f',DCM.Pp.Bs{1}(j,i)));    
        end
    end
    axis ij
    axis off
    

case {lower('Sin(Data) - Region 1')}
    plot_trials(DCM,1);
   
case {lower('Sin(Data) - Region 2')}
    plot_trials(DCM,2);

case {lower('Sin(Data) - Region 3')}
    if nr>2
        plot_trials(DCM,3);
    else
        disp('There is no region 3');
    end
case {lower('Sin(Data) - Region 4')}
    if nr>3
        plot_trials(DCM,4);
    else
        disp('There is no region 4');
    end

end


function [] = plot_trials(DCM,j)

    % Only show a maximum of max_nt trials
    max_nt=16;
    nt     = length(DCM.xY.y);     
    nt=min(nt,max_nt);
    
    rnt=ceil(sqrt(nt));
    for i=1:nt,
        subplot(rnt,rnt,i);
        plot(DCM.xY.pst,sin(DCM.xY.y{i}(:,j)));
        hold on
        plot(DCM.xY.pst,sin(DCM.y{i}(:,j)),'r');
        title(sprintf('Trial %d: %s',i,DCM.xY.code{i}));
    end
drawnow
