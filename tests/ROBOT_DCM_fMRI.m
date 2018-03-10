function E = ROBOT_DCM_fMRI
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

% $Id: ROBOT_DCM_fMRI.m 7279 2018-03-10 21:22:44Z karl $

% tests of spatial models: 'ECD', 'LFP' or 'IMG'
%==========================================================================
try
    cd(fullfile(spm('Dir'),'tests','data','fMRI'));
catch
    try
        cd('C:\Users\karl\Documents\SPM\DCM tests');
        
    catch
        cd('C:\home\spm\DCM\DCM tests');
    end
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

OPT   = {'standard','two_state','stochastic'};
for i = 1:length(OPT)
    try
        % report
        %------------------------------------------------------------------
        fprintf('\nTesting %s:\n\n',OPT{i})
        
        TCM = DCM;
        eval(['TCM.options.' OPT{i} ' = 1;']);
        TCM = spm_dcm_estimate(TCM);
        
        % diagnositics
        %------------------------------------------------------------------
        spm_dcm_fmri_check(TCM);
        set(gcf,'Name',['DCM Diagnositics: Option: ' OPT{i}]);
        
        % print graphics
        %------------------------------------------------------------------
        spm_demo_print
        fprintf('\n\n     --------***--------   \n\n')
        
        % print graphics
        %------------------------------------------------------------------
        F(i) = TCM.F;
        
        if TCM.options.two_state
            a      = exp(spm_vec(TCM.Ep.A))/8;
            A(:,i) = a - 2*diag(diag(a));
            B(:,i) = exp(spm_vec(TCM.Ep.B))/8;
        else
            a      = TCM.Ep.A;
            a      = a - diag(diag(a)) - diag(exp(diag(a))/2);
            A(:,i) = spm_vec(a);
            B(:,i) = spm_vec(TCM.Ep.B);
        end
        
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
title('log-evidence','FontSize',16)
set(gca,'XTickLabel',OPT)
axis square

subplot(2,2,2)
imagesc(A)
title('MAP estimates (A)','FontSize',16)
set(gca,'XTickLabel',OPT)
axis square

subplot(2,2,3)
bar(A)
title('MAP estimates (A)','FontSize',16)
axis square
legend(OPT)

subplot(2,2,4)
bar(B(find(B(:,1)),:))
title('MAP estimates (B)','FontSize',16)
axis square


% print graphics and save
%------------------------------------------------------------------
spm_demo_print

save ROBOT


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
