function E = DCM_fMRI_ROBOT
% test routine to check current implementations of DCM for fMRI
%==========================================================================
%
% Options
%--------------------------------------------------------------------------
% DCM.options.two_state              % two regional populations (E and I)
% DCM.options.stochastic             % fluctuations on hidden states
% DCM.options.nonlinear              % interactions among hidden states
% DCM.options.nograph                % graphical display
% DCM.options.centre                 % mean-centre inputs
% DCM.options.P                      % starting estimates for parameters
% DCM.options.hidden                 % indices of hidden regions


% tests of spatial models: 'ECD', 'LFP' or 'IMG'
%==========================================================================
try
    cd('C:\home\spm\DCM\DCM fMRI')
catch
    cd('C:\Users\karl\Documents\SPM\DCM fMRI')
end
close all
delete(get(0,'Children'))
if exist('DEMO.ps','file')
    delete('DEMO.ps')
end
clc
E = {};


% spatial models
%==========================================================================
load DCM_attention
DCM.options.two_state  = 0;        % two regional populations (E and I)
DCM.options.stochastic = 0;        % fluctuations on hidden states
DCM.options.nonlinear  = 0;        % interactions among hidden states
DCM.options.centre     = 0;        % mean-centre inputs
DCM.options.hidden     = [];       % indices of hidden regions

OPT = {'two_state','stochastic','nonlinear','centre'};
for i = 1:length(OPT)
    try
        % report
        %------------------------------------------------------------------
        fprintf('\ntesting %s\n',OPT{i})
        
        eval(['DCM.options.' OPT{i} ' = 1;']);
        TCM = spm_dcm_estimate(DCM);
        
        spm_figure('GetWin',['fMRI model: ' OPT{i}]);
        spm_dcm_fmri_check(TCM);
        
        % print graphics
        %------------------------------------------------------------------
        spm_demo_print
        
        % print graphics
        %------------------------------------------------------------------
        F(i) = TCM.F;        
        
    catch
        
        % errors
        %------------------------------------------------------------------
        E{end + 1} = lasterror;
        
    end
end

fprintf('\n\n     --------***--------   \n\n')
    
spm_figure('GetWin','Model comparison ERP');

subplot(2,2,1)
bar(F - min(F))
ylabel('log-evidence','FontSize',16)
set(gca,'XTickLabel',str)
axis square


return


% Tests of neuronal models: 'ERP', 'SEP', 'CMC', 'NMM' or 'MFM'
%==========================================================================
load DCM_MMN
DCM.options.spatial  = 'ECD';
DCM.options.analysis = 'ERP';
DCM.name             = 'DCM_MMN';

% neural models
%--------------------------------------------------------------------------
model = {'ERP','SEP','LFP','CMC','CMM','NMM','MFM'};
for i = 1:length(model)
    
    % report
    %----------------------------------------------------------------------
    fprintf('\nChecking neural models (ERP) %s\n',model{i})
    
    try
        
        % invert model
        %------------------------------------------------------------------
        DCM.options.model  = model{i};
        if isfield(DCM,'M')
            DCM  = rmfield(DCM,'M');
        end
        DCM  = spm_dcm_erp(DCM);
        
        spm_figure('GetWin',['ERP model: ' model{i}]);
        spm_dcm_erp_results(DCM,'ERPs (mode)',gcf);
        
        
        % evidence and cod
        %------------------------------------------------------------------
        F(i)   = DCM.F;
        str{i} = model{i};
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

% compare ERP models
%--------------------------------------------------------------------------
spm_figure('GetWin','Model comparison ERP');

subplot(2,2,1)
bar(F - min(F))
ylabel('log-evidence','FontSize',16)
set(gca,'XTickLabel',str)
axis square

subplot(2,2,2)
bar(R)
ylabel('Residual SSQ','FontSize',16)
set(gca,'XTickLabel',str)
legend({'condition 1','condition 2'})
axis square

spm_demo_print


% test of steady state (CSD) models (and LFP spatial model)
%==========================================================================
load DCM_CSD
DCM.options.spatial  = 'LFP';
DCM.options.analysis = 'CSD';
DCM.name             = 'DCM_CSD';

% neural models
%--------------------------------------------------------------------------
clear F R
model = {'ERP','SEP','LFP','CMC','CMM','NMM','MFM'};
for i = 1:length(model)
    
    % report
    %----------------------------------------------------------------------
    fprintf('\nChecking neural models (CSD) %s\n',model{i})
    
    try
        
        % invert model
        %------------------------------------------------------------------
        DCM.options.model  = model{i};
        if isfield(DCM,'M')
            DCM  = rmfield(DCM,'M');
        end
        DCM    = spm_dcm_csd(DCM);
        spm_figure('GetWin',['CSD model: ' model{i}]);
        spm_dcm_csd_results(DCM,'Cross-spectra (channels)',gcf)
        
        % print graphics
        %------------------------------------------------------------------
        spm_demo_print
        
        % evidence and cod
        %------------------------------------------------------------------
        F(i)   = DCM.F;
        str{i} = model{i};
        
    catch
        
        % errors
        %------------------------------------------------------------------
        E{end + 1} = lasterror;
        
    end
    
    fprintf('\n\n     --------***--------   \n\n')
    
end



% compare CSD models
%--------------------------------------------------------------------------
spm_figure('GetWin','Model comparison CSD');

subplot(2,2,1)
bar(F - min(F))
ylabel('log-evidence','FontSize',16)
set(gca,'XTickLabel',str)
axis square

spm_demo_print


% test of induced response models
%==========================================================================
load DCM_FACES

DCM.options.spatial  = 'ECD';
DCM.options.analysis = 'IND';
DCM.name             = 'DCM_FACES';

fprintf('\nChecking spm_dcm_ind\n')

try
    
    if isfield(DCM,'M')
        DCM  = rmfield(DCM,'M');
    end
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
DCM.options.analysis = 'TFM';
DCM.name             = 'DCM_TFM';

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
    try
        disp(E{i}.stack(end - 1))
        disp(E{i}.stack(1))
    end
    disp('------------------------------------------------')
end



return

% Print sub-function
%==========================================================================
function spm_demo_print

% print graphics
%--------------------------------------------------------------------------
drawnow

H     = sort(get(0,'Children'));
for j = 1:length(H);
    
    figure(H(j))
    axes('position',[.05 .98 .9 .02]);
    text(0,0.5,get(gcf,'Name'),'Fontsize',10,'Fontweight','Bold')
    axis off
    
    spm_print('DEMO.ps',gcf)
    
end
delete(H)

return
