function E = DCM_ROBOT
% test routine to check current impleentations of DCM for electophsyiology
%==========================================================================
%   options.analysis     - 'ERP','CSD', 'IND' or 'TFM
%   options.model        - 'ERP','SEP','CMC','LFP','NNM' or 'MFM'
%   options.spatial      - 'ECD','LFP' or 'IMG'


% tests of spatial models: 'ECD', 'LFP' or 'IMG'
%==========================================================================
cd('C:\home\spm\DCM\DCM tests')
clear all
close all
delete(get(0,'Children'))
if exist('DEMO.ps','file')
    delete('DEMO.ps')
end
clc
E = {};


% spatial models
%==========================================================================
load DCM_MMN
DCM.options.analysis = 'ERP';
DCM.options.model    = 'ERP';


model = {'ECD','IMG'};
for i = 1:length(model)
    
    % report
    %----------------------------------------------------------------------
    fprintf('\nChecking spatial models %s\n',model{i})
    
    try
        % invert model
        %------------------------------------------------------------------
        DCM.options.spatial = model{i};
        DCM  = rmfield(DCM,'M');
        DCM  = spm_dcm_erp(DCM);
        
        spm_figure('GetWin',['Standard ERP model: ' model{i}]);
        spm_dcm_erp_results(DCM,'ERPs (mode)',gcf);
        
        % print graphics
        %------------------------------------------------------------------
        spm_demo_print
        
    catch
        
        % errors
        %------------------------------------------------------------------
        E{end + 1} = lasterror;
        
    end
    
    fprintf('\n\n     --------***--------   \n\n')
    
end


% Tests of neuronal models: 'ERP', 'SEP', 'CMC', 'NMM' or 'MFM'
%==========================================================================
load DCM_MMN
DCM.options.spatial  = 'ECD';
DCM.options.analysis = 'ERP';

% neural models
%--------------------------------------------------------------------------
model = {'ERP','SEP','LFP','CMC','CMM','NMM','MFM'};
for i = 1:length(model)
    
    % report
    %----------------------------------------------------------------------
    fprintf('\nChecking neural models %s\n',model{i})
    
    try
        
        % invert model
        %------------------------------------------------------------------ 
        DCM.options.model  = model{i};
        DCM  = rmfield(DCM,'M');
        DCM  = spm_dcm_erp(DCM);
        
        spm_figure('GetWin',['Standard ECD model: ' model{i}]);
        spm_dcm_erp_results(DCM,'ERPs (mode)',gcf);
        
        
        % print graphics
        %------------------------------------------------------------------
        F(i)  = DCM.F;
        for j = 1:length(DCM.R)
            R(i,j) = std(spm_vec(DCM.R{j}));
        end
        
        % print graphics
        %------------------------------------------------------------------
        spm_demo_print
    catch
        
        % errors
        %------------------------------------------------------------------
        E{end + 1} = lasterror;
        
    end
    
    fprintf('\n\n     --------***--------   \n\n')
end

% compare models
%--------------------------------------------------------------------------
spm_figure('GetWin','Model comparison');

subplot(2,2,1)
bar(F - min(F))
ylabel('log-evidence','FontSize',16)
set(gca,'XTickLabel',model)
axis square

subplot(2,2,2)
bar(R)
ylabel('Residual SSQ','FontSize',16)
set(gca,'XTickLabel',model)
legend({'condition 1','condition 2'})
axis square

spm_demo_print

% test of steady state (CSD) models (and LFP spatial model)
%==========================================================================
load DCM_CSD
DCM.options.model    = 'CMC';
DCM.options.spatial  = 'LFP';
DCM.options.analysis = 'CSD';
fprintf('\nChecking spm_dcm_csd\n')

try
    
    DCM  = rmfield(DCM,'M');
    DCM  = spm_dcm_csd(DCM);
    spm_figure('GetWin','Cross-spectral density under Steady state');
    spm_dcm_csd_results(DCM,'Cross-spectra (channels)',gcf)
    
    % print graphics
    %----------------------------------------------------------------------
    spm_demo_print
catch
    
    % errors
    %----------------------------------------------------------------------
    E{end + 1} = lasterror;
    
end

fprintf('\n\n     --------***--------   \n\n')

% test of induced response models
%==========================================================================
load DCM_FACES

DCM.options.spatial  = 'ECD';
DCM.options.analysis = 'IND';
fprintf('\nChecking spm_dcm_ind\n')
try
    
    DCM  = rmfield(DCM,'M');
    DCM  = spm_dcm_ind(DCM);
    spm_figure('GetWin','Induced responses - 2 conditions');
    spm_dcm_ind_results(DCM,'Time-modes',gcf);
    
    % print graphics
    %----------------------------------------------------------------------
    spm_demo_print
    
catch
    
    % errors
    %----------------------------------------------------------------------
    E{end + 1} = lasterror;
    
end

fprintf('\n\n     --------***--------   \n\n')

% test of time-frequency models
%==========================================================================
load DCM_TFM

DCM.options.spatial  = 'ECD';
DCM.options.analysis = 'IND';
fprintf('\nChecking spm_dcm_tfm\n')

try
    
    DCM  = rmfield(DCM,'M');
    DCM  = spm_dcm_tfm(DCM);
    
    spm_figure('GetWin','induced and evoked responses');
    spm_dcm_tfm_results(DCM,'induced and evoked responses',gcf);
    spm_figure('GetWin','induced and evoked predictions');
    spm_dcm_tfm_results(DCM,'induced and evoked predictions',gcf);
    
    % print graphics
    %----------------------------------------------------------------------
    spm_demo_print
    
catch
    
    % errors
    %----------------------------------------------------------------------
    E{end + 1} = lasterror;
    
end

fprintf('\n\n     --------***--------   \n\n')

% Show failed routines
%--------------------------------------------------------------------------
for i = 1:length(E)
    disp(E{i}.message)
    disp(E{i}.stack(end - 1))
    disp(E{i}.stack(1))
    disp('------------------------------------------------')
end



return

% Print subfunction
%==========================================================================
function spm_demo_print

% print graphics
%--------------------------------------------------------------------------
H     = sort(get(0,'Children'));
for j = 1:length(H);
    
    figure(H(j))
    axes('position',[.05 .98 .9 .02]);
    text(0,0.5,get(gcf,'Name'),'Fontsize',10,'Fontweight','Bold')
    axis off
    
    print(gcf,'-dpsc','-append','-r0','DEMO.ps')
    
end
delete(H)

return
