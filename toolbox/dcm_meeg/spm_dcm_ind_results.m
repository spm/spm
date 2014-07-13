function [DCM] = spm_dcm_ind_results(DCM,Action,fig)
% Results for induced Dynamic Causal Modelling (DCM)
% FORMAT [DCM] = spm_dcm_ind_results(DCM,Action)
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
%     'Input (C x u)'
%     'Dipoles'
%     'Save results as img'           
%__________________________________________________________________________
%
% DCM is a causal modelling procedure for dynamical systems in which
% causality is inherent in the differential equations that specify the
% model.  The basic idea is to treat the system of interest, in this case
% the brain, as an input-state-output system.  By perturbing the system
% with known inputs, measured responses are used to estimate various
% parameters that govern the evolution of brain states.  Although there are
% no restrictions on the parameterisation of the model, a bilinear
% approximation affords a simple re-parameterisation in terms of effective
% connectivity.  This effective connectivity can be latent or intrinsic or,
% through bilinear terms, model input-dependent changes in effective
% connectivity. Parameter estimation proceeds using fairly standard
% approaches to system identification that rest upon Bayesian inference.
%__________________________________________________________________________
% Copyright (C) 2007-2013 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_ind_results.m 6101 2014-07-13 21:34:34Z karl $
 
 
% set up
%--------------------------------------------------------------------------
if nargin < 3
    spm_figure('GetWin','Graphics');
else
    figure(fig);
end
colormap(gray), clf

xY     = DCM.xY;
xU     = DCM.xU;
nt     = length(xY.y);           % Nr of trial types
nr     = size(xY.xf,2);          % Nr of sources
nu     = size(xU.X, 2);          % Nr of experimental effects
nf     = DCM.options.Nmodes;     % Nr of frequency modes explained
ns     = size(xY.y{1},1);        % Nr of time bins
pst    = xY.pst;                 % peri-stmulus time
Hz     = xY.Hz;                  % frequencies
xY.U   = xY.U(:,1:nf);           % remove unmodelled frequency modes
 
    
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
            TF{i,j} = TF{i,j} + xY.xf{i,j}(:,1:nf)*xY.U';
        end
    end
    
    % loop over trials, sources (predicted and observed)
    %----------------------------------------------------------------------
    for i = 1:nt
        for j = 1:nr
           
            subplot(2*nt,nr,(i - 1)*nr + j)
            imagesc(pst,Hz,TF{i,j}')
            axis xy
            xlabel('pst (ms)')
            ylabel('frequency')
            title(sprintf('trial %i: %s ',i,DCM.Sname{j}));
            
            
            % write images if requested
            %--------------------------------------------------------------
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
                V.fname = [sprintf('%s_TFR%d%d',DCM.name(1:end-4),i,j),spm_file_ext];
                spm_write_vol(V, TF{i,j}');
            end
        end
    end
 
case{lower('Frequency modes')}
    
    % spm_dcm_ind_results(DCM,'Frequency modes')
    %----------------------------------------------------------------------
    plot(xY.Hz,xY.U)
    xlabel('Frequnecy (Hz)')
    xlabel('modes')
    title('Frequency modes modelled at each source','FontSize',16)
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
    
    % reconstitute time-frequency and get principal mode over channels
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
    cmax  = zeros(nt, nr);
    for i = 1:nt
        for j = 1:nr
            subplot(nt*2,nr,(i - 1)*2*nr + j)
            imagesc(pst,Hz,TF{i,j}' + RF{i,j}')
            axis xy
            xlabel('pst (ms)')
            ylabel('frequency')
            title({sprintf('trial %i: %s ',i,DCM.Sname{j});
                  'observed (adjusted)'})
              
            clim = caxis;
            cmax(i, j) = max(clim);  
 
            subplot(nt*2,nr,(i - 1)*2*nr + nr + j)
            imagesc(pst,Hz,TF{i,j}')
            axis xy
            xlabel('pst (ms)')
            ylabel('frequency')
            title({sprintf('trial %i',i); 'predicted'})
        end
    end
 
    cmax  = mean(cmax, 1);
    for i = 1:nt
        for j = 1:nr
            subplot(nt*2,nr,(i - 1)*2*nr + j)
            caxis(cmax(j)*[-1 1]);
 
            subplot(nt*2,nr,(i - 1)*2*nr + nr + j)
            caxis(cmax(j)*[-1 1]);
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
                V.fname = [sprintf('%s_A%d%d',DCM.name(1:end-4),i,j),spm_file_ext];
                spm_write_vol(V,A);
            end
 
        end
    end
    
    axes('position', [0.4, 0.95, 0.2, 0.01]);
    axis off;
    title('endogenous coupling (A)','FontSize',16)
    colormap(jet);
     
case{lower('Coupling (B - Hz)')}
    
    % get experimental effect (if any)
    %----------------------------------------------------------------------
    if     nu == 0
        return
    elseif nu == 1;
        k = 1;
    else
        [k, ok] = listdlg('PromptString', 'which effect', 'Name', 'please select',...
            'SelectionMode','single', 'ListString', DCM.xU.name);
        if ~ok, return; end
    end
    
    % reconstitute time-frequency coupling
    %----------------------------------------------------------------------
    for i = 1:nr
        for j = 1:nr
            subplot(nr,nr,j + nr*(i - 1))
            ii = (1:nf)*nr - nr + i;
            jj = (1:nf)*nr - nr + j; 
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
               V.fname = [sprintf('%s_B%d%d',DCM.name(1:end-4),i,j),spm_file_ext];
               spm_write_vol(V,B);
            end
                
        end
    end
    
    axes('position', [0.4, 0.95, 0.2, 0.01]);
    axis off;
    title({'changes in coupling (B)';xU.name{k}});
    colormap(jet);
    
case{lower('Coupling (A - modes)')}
    
        
    % images
    %----------------------------------------------------------------------
    subplot(3,1,1)
    imagesc(DCM.Ep.A)

    set(gca,'YTick',(1:nr),'YTickLabel',DCM.Sname,'FontSize',8)
    set(gca,'XTick',[])
    xlabel('from','FontSize',10)
    ylabel('to','FontSize',10)
    title('Coupling (Hz)','FontSize',16)
    axis square
    colorbar
 
 
    % PPM
    %----------------------------------------------------------------------
    subplot(3,1,2)
    imagesc(DCM.Pp.A)
    set(gca,'YTick',(1:nr),'YTickLabel',DCM.Sname,'FontSize',8)
    set(gca,'XTick',[])
    title('Conditional probabilities','FontSize',16)
    axis square
    colorbar   
    
    % Guide
    %----------------------------------------------------------------------
    subplot(3,2,5)
    image(48*(kron(eye(nf,nf),ones(nr,nr)) - speye(nr*nf,nr*nf)))
    title('Within frequency (linear)','FontSize',12)
    axis square
    
    subplot(3,2,6)
    image(48*(kron(1 - eye(nf,nf),ones(nr,nr))))
    title('Between frequency (non-linear)','FontSize',12)
    axis square
 
 
case{lower('Coupling (B - modes)')}
    
    % spm_dcm_erp_results(DCM,'coupling (B)');
    %----------------------------------------------------------------------
    for i = 1:min(nu,4)
        
        % images
        %------------------------------------------------------------------
        subplot(4,2,2*(i - 1) + 1)
        imagesc(DCM.Ep.B{i})
        set(gca,'YTick',(1:nr),'YTickLabel',DCM.Sname,'FontSize',8)
        set(gca,'XTick',[])
        xlabel('from')
        ylabel('to')
        title([DCM.xU.name{i} ' (Hz)'],'FontSize',16)
        axis square
        colorbar
 
        
        % PPM
        %------------------------------------------------------------------
        subplot(4,2,2*(i - 1) + 2)
        imagesc(DCM.Pp.B{i})
        set(gca,'YTick',(1:nr),'YTickLabel',DCM.Sname,'FontSize',8)
        set(gca,'XTick',[])
        title('Conditional probabilities','FontSize',16)
        axis square
        colorbar
 
    end
 
 
case{lower('Input (C - Hz)')}
   
    
    % reconstitute time-frequency and get principal mode over channels
    %----------------------------------------------------------------------
    nu    = size(DCM.Ep.C,2);
    for k = 1:nu
        subplot(nu,1,k)
        for i = 1:nr
            j = [1:nf]*nr - nr + i;
            UF(:,i) = xY.U*DCM.Ep.C(j,k);
        end
        plot(Hz,UF)
        xlabel('Frequency (Hz)')
        title(sprintf('frequency response to input %i',k),'FontSize',16)
        axis square, grid on
    end
    legend(DCM.Sname)
    
    
case{lower('Input (u - ms)')}
    
    % get input
    % ---------------------------------------------------------------------
    U    = spm_erp_u((pst - pst(1))/1000,DCM.Ep,DCM.M);
    
    subplot(1,1,1)
    plot(pst,U)
    xlabel('time (ms)')
    title('input','FontSize',16)
    axis square tight, grid on
    
case{lower('Input (C x u)')}
    
    % get input
    % ---------------------------------------------------------------------
    U    = spm_erp_u((pst - pst(1))/1000,DCM.Ep,DCM.M);
    U    = DCM.Ep.C*U';
    
    for k = 1:nr
        
        % reconstitute time-frequency input
        %------------------------------------------------------------------
        j  = [1:nf]*nr - nr + k;
        UF = xY.U*U(j,:);
        
        subplot(nr,1,k)
        imagesc(pst,Hz,UF)
        xlabel('PST (ms)')
        ylabel('Frequency (Hz)')
        title(sprintf('input to %s',DCM.Sname{k}),'FontSize',16)
        axis xy
        
    end
    
    
case{lower('Dipoles')}
    
    sdip.n_seeds = 1;
    sdip.n_dip  = nr;
    sdip.Mtb    = 1;
    sdip.j{1}   = zeros(3*nr, 1);
    sdip.loc{1} = full(DCM.M.dipfit.Lpos);
    spm_eeg_inv_ecd_DrawDip('Init', sdip)
        
 
case{lower('Save results as img')}
 
    fprintf('Saving the Time-frequency representation at sources\n');
    DCM.saveInd = 'TFR';
    spm_dcm_ind_results(DCM,'Wavelet');
 
    fprintf('Saving the coupling matrix A\n');
    DCM.saveInd ='Amatrix';
    spm_dcm_ind_results(DCM,'Coupling (A - Hz)');
 
    fprintf('Saving the coupling matrix B\n');
    DCM.saveInd = 'Bmatrix';
    spm_dcm_ind_results(DCM,'Coupling (B - Hz)');
    DCM = rmfield(DCM,'saveInd');
 
end
drawnow
