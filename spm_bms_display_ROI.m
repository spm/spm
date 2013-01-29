function spm_bms_display_ROI (BMS,mask,method)
% display results from BMS in a region of interest (ROI)
% FORMAT spm_bms_display_ROI (BMS,mask,method)
%
% Input:
% BMS    - BMS.mat file 
% mask   - region of interest image
% method - inference method (FFX or RFX)
%__________________________________________________________________________
% Copyright (C) 2009-2011 Wellcome Trust Centre for Neuroimaging

% Maria Joao Rosa
% $Id: spm_bms_display_ROI.m 5219 2013-01-29 17:07:07Z spm $

% Find graphics window
% -------------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');

% Input
% -------------------------------------------------------------------------
if nargin<1
    BMS    = spm_select(1,'^BMS.mat$','select BMS.mat files');
end
if nargin<2
    mask   = spm_select(1,'image','select ROI file');
end
if nargin<3
    method = spm_input('Inference method',+1,'b','FFX|RFX',['ffx';'rfx']);
end

mask_image = spm_vol(mask);         % Mask image Vol

% Nb. of subjects and models
% -------------------------------------------------------------------------
if isfield(BMS.map,method)
   switch method
       case 'ffx'
            data   = BMS.map.ffx.data;
       case 'rfx'
            data   = BMS.map.rfx.data;
   end
else
  msgbox(sprintf('Error: no %s analysis in current BMS.mat!',method))
  return                  
end

nsubjs    = size(data,2);
nmodels   = size(data{1}(1).mod_map,1);
nsess     = size(data{1},2);

% Sort out log-evidence images dimensions
% -------------------------------------------------------------------------
Vol_models(1,1) = spm_vol(data{1}.mod_map(1));

first_vol       = Vol_models(1,1);
M               = first_vol{1}.mat;
DIM             = first_vol{1}.dim(1:3)'; 

xdim            = DIM(1); 
ydim            = DIM(2); 
zdim            = DIM(3);
[xords,yords]   = ndgrid(1:xdim,1:ydim);
xords           = xords(:)';  
yords           = yords(:)';
I               = 1:xdim*ydim;
zords_init      = ones(1,xdim*ydim);

% Loop through data
% -------------------------------------------------------------------------
for i = 1:nmodels,
    for s = 1:nsubjs,
        for se = 1:nsess,
             Vol_models(s,i,se) = spm_vol(data{s}(se).mod_map(i));
        end
    end
end

log_ev_roi = zeros(nsubjs,nmodels);
nvox_total = 0;

% Loop through image slices
% -------------------------------------------------------------------------
for z = 1:zdim,
    
    j = NaN(xdim,ydim);                         % Init. image values
 
    zords   = z*zords_init;                     % Slice z
    xyz     = [xords(I); yords(I); zords(I)];   % Slice coordinates
    nVox    = size(xyz,2);                      % Nb. of voxels per slice
    
    % Voxels inside mask
    mask_xyz  = mask_image.mat\M*[xyz(:,1:nVox);ones(1,nVox)];
    gamma     = spm_get_data(mask_image,mask_xyz);
    b         = find(gamma>0.5);            % Voxels in the mask
    
    z_models        = NaN(nsubjs,nmodels,nVox);       % Data 
    z_models(1,1,:) = spm_get_data(first_vol{1},xyz); % Data: all subs/mods  
    non_nan         = find(~isnan(z_models(1,1,:)));  % Voxels ~NaN

    % Find voxels ~NaN and sum sessions
    % ---------------------------------------------------------------------
    for s = 1:nsubjs,
        for k = 1:nmodels,
                sum_tmp_data    = [];
            for ns = 1:nsess,
                tmp_data        = Vol_models(s,k,ns);
                sum_tmp_data    = [sum_tmp_data; spm_get_data(tmp_data{1},xyz)];
            end
                z_models(s,k,:) = sum(sum_tmp_data,1);
                non_nani        = find(~isnan(z_models(s,k,:)));
                non_nan         = intersect(non_nan,non_nani);
        end
    end

    % Voxels to be analysed
    non_nan = intersect(non_nan,b);    
    Nvoxels = length(non_nan);
    
    if Nvoxels > 0
        
        nvox_total = nvox_total + Nvoxels;
        % Do BMS in all voxels of slice z
        for n = 1:Nvoxels,
            log_ev_tmp = z_models(:,:,non_nan(n));
            % Group BMS
            log_ev_roi = log_ev_roi + log_ev_tmp;

        end
    end
end % Loop over slices

% Method
% -------------------------------------------------------------------------
switch method
    
    % Fixed Effects
    % ---------------------------------------------------------------------
    case 'ffx',
        
        F    = sum(log_ev_roi) - min(log_ev_roi);
        i    = F < (max(F) - 32);
        P    = F;
        P(i) = max(F) - 32;
        P    = P - min(P);
        P    = exp(P);
        P    = P/sum(P);
        
        if  isfield(BMS.map,'ffx')
                      
            % Bar plot
            figure(Fgraph);
            spm_results_ui('Clear',Fgraph);
            
            hvox   = axes('Position',[0.25 0.15 0.5 0.25],'Parent',...
                Fgraph,'Visible','off');
            
            bar(1:nmodels,P)
            set(gca,'XTick',1:nmodels)
            set(gca,'XTickLabel',1:nmodels)
            set(gca,'YLim',[0 1])
            ylabel('Posterior Model Probability','Fontsize',12)
            xlabel('Models','Fontsize',12)
            title({'Fixed-effects BMS';''},...
                'Fontsize',12);
            axis square
            grid on
            
            return
            
        else
            
            msgbox('Error: no FFX analysis in current BMS.mat!')
            return
            
        end
        
        
    % Random Effects
    % ---------------------------------------------------------------------
    case 'rfx',
        
        nsamps           = 1e3;
        [alpha,exp_r,xp] = spm_BMS(log_ev_roi,nsamps,0,0,1);
        
        if  isfield(BMS.map,'rfx')
        
            % Bar plots       
            figure(Fgraph);
            spm_results_ui('Clear',Fgraph); 
                
            hvox   = axes('Position',[0.55 0.18 0.30 0.20],'Parent',...
                    Fgraph,'Visible','off');
        
            bar(1:nmodels,xp)
            set(gca,'XTick',1:nmodels)
            set(gca,'XTickLabel',1:nmodels)
            set(gca,'YLim',[0 1])
            ylabel('Exceedance Probability','Fontsize',12)
            xlabel('Models','Fontsize',12)
            title({'Random-effects BMS';''},'Fontsize',12)
            axis square
            grid on
            
            hvox   = axes('Position',[0.16 0.18 0.30 0.20],'Parent',...
                    Fgraph,'Visible','off'); 

            bar(1:nmodels,exp_r)
            set(gca,'XTick',1:nmodels)
            set(gca,'XTickLabel',1:nmodels)
            set(gca,'YLim',[0 1])
            ylabel('Expected Posterior Probability','Fontsize',12)
            xlabel('Models','Fontsize',12)
            title({'Random-effects BMS';''},'Fontsize',12)
            axis square
            grid on

            return
        
        else
                
            msgbox('Error: no RFX analysis in current BMS.mat!')
            return

        end
        
end