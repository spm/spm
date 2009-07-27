function spm_bms_display_vox(BMS,xyz)
% display results from BMS Maps at current voxel
% FORMAT spm_bms_display_vox(xyz)
%
% Input:
% xyz - voxel coordinates [1,3] (voxel)
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Maria Joao Rosa
% $Id: spm_bms_display_vox.m 3288 2009-07-27 09:23:54Z maria $

% Find graphics window
% -------------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');

% Inference method to plot
% -------------------------------------------------------------------------
method = spm_input('Inference method',+1,'b','FFX|RFX',['FFX';'RFX']);

switch method
    
    % Fixed effects
    % ---------------------------------------------------------------------
    case 'FFX'
        
            if  isfield(BMS.map,'ffx')
                
                nmodels = size(BMS.map.ffx.ppm,2);
                models  = [];
                ppm_vox = zeros(nmodels,1);
        
                % Get values
                for i = 1:nmodels,
                    tmp_ppm_vox   = spm_vol(BMS.map.ffx.ppm{i});
                    ppm_vox(i,:)  = spm_get_data(tmp_ppm_vox,xyz);
                    models        = [models; sprintf('model %d',i)];
                end
         
                log_gbf = log(ppm_vox./(1-ppm_vox));
                
                % Bar plot
                figure(Fgraph);
                spm_results_ui('Clear',Fgraph);
        
                hvox   = axes('Position',[0.25 0.15 0.5 0.25],'Parent',...
                Fgraph,'Visible','off');

                bar(1:nmodels,ppm_vox)
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

        

     
    % Random effects
    % ---------------------------------------------------------------------
    case 'RFX'
        
        if  isfield(BMS.map,'rfx')
            
            nmodels    = size(BMS.map.rfx.ppm,2);      
            models     = [];
            exp_r_vox  = zeros(nmodels,1);
            epm_exists   = isfield(BMS.map.rfx,'epm');
            
            if  epm_exists,   xp_vox    = zeros(nmodels,1); end
        
            % Get values
            for i = 1:nmodels,
                tmp_exp_r_vox   = spm_vol(BMS.map.rfx.ppm{i});
                exp_r_vox(i,:)  = spm_get_data(tmp_exp_r_vox,xyz);
                if epm_exists
                tmp_xp_vox      = spm_vol(BMS.map.rfx.epm{i});
                xp_vox(i,:)     = spm_get_data(tmp_xp_vox,xyz);
                end
                models          = [models; sprintf('model %d',i)];
            end
        
            % Bar plots       
            figure(Fgraph);
            spm_results_ui('Clear',Fgraph); 
        
            if epm_exists,
                
            hvox   = axes('Position',[0.55 0.18 0.30 0.20],'Parent',...
                    Fgraph,'Visible','off');
        
            bar(1:nmodels,xp_vox)
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
                
            else
                
            hvox   = axes('Position',[0.25 0.15 0.5 0.25],'Parent',...
                Fgraph,'Visible','off');
            
            end

            bar(1:nmodels,exp_r_vox)
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

end