function [DCM] = spm_dcm_ind_results(DCM,Action)
% Results for induced Dynamic Causal Modeling (DCM)
% FORMAT [DCM] = spm_dcm_ind_results(DCM,Action);
% Action:
%     'Frequency modes'
%     'Time-modes'
%     'Time-frequency'
%     'Coupling (A - Hz)'
%     'Coupling (B - Hz)'
%     'Coupling (A - modes)'
%     'Coupling (B - modes)'
%     'Input (C - Hz)'
%     'Input (u - ms)'
%     'Dipoles'
%     'Saveimg'           
%__________________________________________________________________________
%
% DCM is a causal modelling procedure for dynamical systems in which
% causality is inherent in the differential equations that specify the model.
% The basic idea is to treat the system of interest, in this case the brain,
% as an input-state-output system.  By perturbing the system with known
% inputs, measured responses are used to estimate various parameters that
% govern the evolution of brain states.  Although there are no restrictions
% on the parameterisation of the model, a bilinear approximation affords a
% simple re-parameterisation in terms of effective connectivity.  This
% effective connectivity can be latent or intrinsic or, through bilinear
% terms, model input-dependent changes in effective connectivity. Parameter
% estimation proceeds using fairly standard approaches to system
% identification that rest upon Bayesian inference.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_ind_results.m 3362 2009-09-04 14:21:40Z guillaume $


% get figure handle
%--------------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');
colormap(gray)
figure(Fgraph)
clf

xY     = DCM.xY;
nt     = length(xY.y);           % Nr of trial types
nr     = length(DCM.C);          % Nr of sources
nu     = length(DCM.B);          % Nr of experimental effects
nf     = size(xY.U,2);           % Nr of frequency modes
ns     = size(xY.y{1},1);        % Nr of time bins
pst    = xY.pst;                 % peri-stmulus time
Hz     = xY.Hz;                  % frequencies

    
% switch
%--------------------------------------------------------------------------
switch(lower(Action))
    
    
    % display time-frequency data if requested
    %----------------------------------------------------------------------
    case{lower('Wavelet')}
    
    % reconstitute time-frequency and get principle model over channels
    %----------------------------------------------------------------------
    nk    = length(Hz);
    TF    = cell(nt,nr);
    for i = 1:nt
        for j = 1:nr
            TF{i,j} = sparse(ns,nk);
        end
    end
    for i = 1:nt
        for j = 1:nr
            TF{i,j} = TF{i,j} + xY.xf{i,j}*xY.U';
        end
    end
    
    % loop over trials, sources (predicted and observed)
    %----------------------------------------------------------------------
    colormap(jet)
    for i = 1:nt
        for j = 1:nr
           
            subplot(2*nt,nr,(i - 1)*nr + j)
            imagesc(pst,Hz,TF{i,j}')
            axis xy
            xlabel('pst (ms)')
            ylabel('frequency')
            title(sprintf('trial %i: %s ',i,DCM.Sname{j}));
            
            if isfield(DCM,'saveInd')&& strcmp(DCM.saveInd,'TFR')
                V.dt    = [spm_type('float64') spm_platform('bigend')];
                V.mat   = [Hz(2)-Hz(1)  0              0  min(Hz);...
                           0            pst(2)-pst(1)  0  min(pst);...
                           0            0              1  0;...
                           0            0              0  1];
                V.mat(1,4) = V.mat(1,4) - V.mat(1,1);
                V.mat(2,4) = V.mat(2,4) - V.mat(2,2);
                V.pinfo = [1 0 0]';
                V.dim   = [length(Hz) length(pst)  1 ];
                V.fname =sprintf('%s_TFR%d%d.img',DCM.name(1:end-4),i,j);
                spm_write_vol(V, TF{i,j}');
            end
        end
    end

case{lower('Frequency modes')}
    
    % spm_dcm_ind_results(DCM,'Frequency modes')
    %----------------------------------------------------------------------
    plot(DCM.xY.Hz,DCM.xY.U)
    xlabel('Frequnecy (Hz)')
    xlabel('modes')
    title('Frequency modes modelled at each source')
    axis square
    grid on
    
    for i = 1:nf
        str{i} = sprintf('mode %i',i);
    end
    legend(str)
    
case{lower('Time-modes')}
    
    
    % spm_dcm_ind_results(DCM,'Time-modes');
    %----------------------------------------------------------------------  
    for i = 1:nf
        subplot(ceil(nf/2),2,i), hold on
        str   = {};
        for j = 1:nt
            plot(pst,DCM.H{j,i},'LineWidth',2);
            plot(pst,DCM.H{j,i} + DCM.R{j,i},'-.');
            set(gca, 'XLim', [pst(1) pst(end)]);
        end
        hold off
        title({sprintf('mode %i (all regions/trials)',i),'- predicted; -- observed'})
        grid on
        axis square
        xlabel('time (ms)')
        if i == 1
            ylim1 = ylim;
        end

        ylim(max(abs(ylim1))*[-1 1]);
            
    end
    legend(DCM.Sname)

    
case{lower('Time-frequency')}
    
    % reconstitute time-frequency and get principle model over channels
    %----------------------------------------------------------------------
    nk    = length(Hz);
    for i = 1:nt
        for j = 1:nr
            TF{i,j} = sparse(ns,nk);
            RF{i,j} = sparse(ns,nk);
        end
    end
    for i = 1:nt
        for j = 1:nr
            for k = 1:nf
                TF{i,j} = TF{i,j} + DCM.H{i,k}(:,j)*xY.U(:,k)';
                RF{i,j} = RF{i,j} + DCM.R{i,k}(:,j)*xY.U(:,k)';
            end
        end
    end

      
    % loop over trials, sources (predicted and observed)
    %----------------------------------------------------------------------
    colormap(jet)
    cmax = zeros(nt, nr);
    for i = 1:nt
        for j = 1:nr
            subplot(nt*2,nr,(i - 1)*2*nr + j)
            imagesc(pst,Hz,TF{i,j}' + RF{i,j}')
            axis xy
            xlabel('pst (ms)')
            ylabel('frequency')
            title({sprintf('trial %i: %s ',i,DCM.Sname{j});
                  'observed (adjusted for confounds)'})
              
            clim = caxis;  
            cmax(i, j) = max(clim);  

            subplot(nt*2,nr,(i - 1)*2*nr + nr + j)
            imagesc(pst,Hz,TF{i,j}')
            axis xy
            xlabel('pst (ms)')
            ylabel('frequency')
            title({sprintf('trial %i',i);
                  'predicted'})
        end
    end

    cmax = mean(cmax);

    for i = 1:nt
        for j = 1:nr
            subplot(nt*2,nr,(i - 1)*2*nr + j)
            caxis(cmax(j)*[0 1]);

            subplot(nt*2,nr,(i - 1)*2*nr + nr + j)
            caxis(cmax(j)*[0 1]);
        end
    end

case{lower('Coupling (A - Hz)')}
    
    % reconstitute time-frequency coupling
    %----------------------------------------------------------------------
    for i = 1:nr
        for j = 1:nr
            subplot(nr,nr,j + nr*(i - 1))
            ii = [1:nf]*nr - nr + i;
            jj = [1:nf]*nr - nr + j; 
            A  = xY.U*DCM.Ep.A(ii,jj)*xY.U';
            imagesc(Hz,Hz,A)
            caxis(max(abs(caxis))*[-1 1]);                
            axis image
            
            % source names
            %--------------------------------------------------------------
            if i == 1, title({'from'; DCM.Sname{j}}), end
            if j == 1, ylabel({'to';  DCM.Sname{i}}), end
            
            if isfield(DCM,'saveInd')&& strcmp(DCM.saveInd,'Amatrix')
                V.dt    = [spm_type('float64') spm_platform('bigend')];
                V.mat   = [Hz(2)-Hz(1)  0            0  min(Hz);...
                           0            Hz(2)-Hz(1)  0  min(Hz);...
                           0            0            1  0;...
                           0            0            0  1];
                V.mat(1,4) = V.mat(1,4) - V.mat(1,1);
                V.mat(2,4) = V.mat(2,4) - V.mat(2,2);
                V.pinfo = [1 0 0]';
                V.dim   = [length(Hz) length(Hz)  1 ];
                V.fname = sprintf('%s_A%d%d.img',DCM.name(1:end-4),i,j);
                spm_write_vol(V,A);
            end

        end
    end
    
    axes('position', [0.4, 0.95, 0.2, 0.01]);
    axis off;
    title('endogenous coupling (A)')
    colormap(jet);
     
case{lower('Coupling (B - Hz)')}
    
    % get experimental effect (if any)
    %----------------------------------------------------------------------
    if     nu == 0
        return
    elseif nu == 1;
        k = 1;
    else
        k = questdlg('which effect','please select',DCM.xU.name{:},DCM.xU.name{1});
        for i = 1:length(DCM.xU.name)
            b(i) = strcmp(k,DCM.xU.name{i});
        end
        k = find(b);
    end
    
    % reconstitute time-frequency  coupling
    %----------------------------------------------------------------------
    for i = 1:nr
        for j = 1:nr
            subplot(nr,nr,j + nr*(i - 1))
            ii = [1:nf]*nr - nr + i;
            jj = [1:nf]*nr - nr + j; 
            B  = xY.U*DCM.Ep.B{k}(ii,jj)*xY.U';
            
            imagesc(Hz,Hz,B)
            caxis(max(abs(caxis))*[-1 1]);  
            axis image
            
            % source names
            %--------------------------------------------------------------
            if i == 1, title({'from'; DCM.Sname{j}}), end
            if j == 1, ylabel({'to';  DCM.Sname{i}}), end

            
            if isfield(DCM,'saveInd')&& strcmp(DCM.saveInd,'Bmatrix')
               V.dt    = [spm_type('float64') spm_platform('bigend')];
               V.mat   = [Hz(2)-Hz(1)  0            0  min(Hz);...
                          0            Hz(2)-Hz(1)  0  min(Hz);...
                          0            0            1  0;...
                          0            0            0  1];
               V.mat(1,4) = V.mat(1,4) - V.mat(1,1);
               V.mat(2,4) = V.mat(2,4) - V.mat(2,2);
               V.pinfo = [1 0 0]';
               V.dim   = [length(Hz) length(Hz)  1 ];
               V.fname = sprintf('%s_B%d%d.img',DCM.name(1:end-4),i,j);
               spm_write_vol(V,B);
            end
                
        end
    end
    
    axes('position', [0.4, 0.95, 0.2, 0.01]);
    axis off;
    title({'changes in coupling (B)';DCM.xU.name{k}});
    colormap(jet);
    
case{lower('Coupling (A - modes)')}
    
        
    % images
    %----------------------------------------------------------------------
    subplot(3,2,1)
    imagesc(DCM.Ep.A)
    title('Coupling','FontSize',10)
    set(gca,'YTick',[1:nr],'YTickLabel',DCM.Sname,'FontSize',8)
    set(gca,'XTick',[])
    xlabel('from','FontSize',10)
    ylabel('to','FontSize',10)
    axis square

    % table
    %----------------------------------------------------------------------
    subplot(3,2,2)
    text(-1/8,1/2,num2str(DCM.Ep.A,' %-8.2f'),'FontSize',8)
    axis off,axis square


    % PPM
    %----------------------------------------------------------------------
    subplot(3,2,3)
    image(64*DCM.Pp.A)
    set(gca,'YTick',[1:nr],'YTickLabel',DCM.Sname,'FontSize',8)
    set(gca,'XTick',[])
    title('Conditional probabilities')
    axis square

    % table
    %----------------------------------------------------------------------
    subplot(3,2,4)
    text(-1/8,1/2,num2str(DCM.Pp.A,' %-8.2f'),'FontSize',8)
    axis off, axis square
    
    
    % Guide
    %----------------------------------------------------------------------
    subplot(3,2,5)
    image(48*(kron(eye(nf,nf),ones(nr,nr)) - speye(nr*nf,nr*nf)))
    title('Within frequency (linear)')
    axis square
    
    subplot(3,2,6)
    image(48*(kron(1 - eye(nf,nf),ones(nr,nr))))
    title('Between frequency (non-linear)')
    axis square


case{lower('Coupling (B - modes)')}
    
    % spm_dcm_erp_results(DCM,'coupling (B)');
    %----------------------------------------------------------------------
    for i = 1:nu
        
        % images
        %------------------------------------------------------------------
        subplot(4,nu,i)
        imagesc(DCM.Ep.B{i})
        title(DCM.xU.name{i},'FontSize',10)
        set(gca,'YTick',[1:ns],'YTickLabel',DCM.Sname,'FontSize',8)
        set(gca,'XTick',[])
        xlabel('from','FontSize',8)
        ylabel('to','FontSize',8)
        axis square

        % tables
        %------------------------------------------------------------------
        subplot(4,nu,i + nu)
        text(-1/8,1/2,num2str(full(DCM.Ep.B{i}),' %-8.2f'),'FontSize',8)
        axis off
        axis square
        
        % PPM
        %------------------------------------------------------------------
        subplot(4,nu,i + 2*nu)
        image(64*DCM.Pp.B{i})
        set(gca,'YTick',[1:ns],'YTickLabel',DCM.Sname,'FontSize',8)
        set(gca,'XTick',[])
        title('PPM')
        axis square

        % tables
        %------------------------------------------------------------------
        subplot(4,nu,i + 3*nu)
        text(-1/8,1/2,num2str(DCM.Pp.B{i},' %-8.2f'),'FontSize',8)
        axis off
        axis square
        
    end


case{lower('Input (C - Hz)')}
   
    
    % reconstitute time-frequency and get principal mode over channels
    %----------------------------------------------------------------------
    Nu    = size(DCM.Ep.C,2);
    for k = 1:Nu
        subplot(Nu,1,k)
        for i = 1:nr
            j = [1:nf]*nr - nr + i;
            UF(:,i) = xY.U*DCM.Ep.C(j,k);
        end
        plot(Hz,UF)
        xlabel('Frequency (Hz)')
        title(sprintf('frequency response to input %i',k))
        axis square, grid on
    end
    legend(DCM.Sname)
    
    
case{lower('Input (u - ms)')}
    
    % get input
    % ---------------------------------------------------------------------
    U    = spm_ind_u((pst - pst(1))/1000,DCM.Ep,DCM.M);
    
    subplot(1,1,1)
    plot(pst,U)
    xlabel('time (ms)')
    title('input')
    axis square tight, grid on
    
    
case{lower('Dipoles')}
    
    sdip.n_seeds = 1;
    sdip.n_dip  = nr;
    sdip.Mtb    = 1;
    sdip.j{1}   = zeros(3*nr, 1);
    sdip.loc{1} = full(DCM.M.dipfit.Lpos);
    spm_eeg_inv_ecd_DrawDip('Init', sdip)
        

case{lower('Save results as img')}

    fprintf('Saving the Time-frequency representation at sources\n');
    DCM.saveInd='TFR';
    spm_dcm_ind_results(DCM,'Wavelet');

    fprintf('Saving the coupling matrix A\n');
    DCM.saveInd='Amatrix';
    spm_dcm_ind_results(DCM,'Coupling (A - Hz)');

    fprintf('Saving the coupling matrix B\n');
    DCM.saveInd='Bmatrix';
    spm_dcm_ind_results(DCM,'Coupling (B - Hz)');
    DCM=rmfield(DCM,'saveInd');

end
drawnow

