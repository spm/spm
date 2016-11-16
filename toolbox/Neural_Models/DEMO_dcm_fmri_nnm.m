function DCM = DEMO_dcm_fmri_nnm
% This demonstration routine illustrates the dynamic causal modelling of
% fMRI timeseries using neural mass models. We first specify a simple DCM
% for the attentional dataset. Following inversion, the posterior densities
% are used to characterise the haemodynamic correlates of induced
% responses.
%
% This experiment involved attention to visual motion. We then use Bayesian
% model reduction to ask whether attention was mediated through the
% modulation of deep or superficial pyramidal cells in the visual motion
% sensitive area (V5 or MST). in this setup, there are three regions and
% three inputs (visual stimulation, visual motion and attention  to
% motion). We treat the latter two inputs as modulatory; namely, increasing
% extrinsic or intrinsic connectivity in particular parts of the network.
% Intrinsic connectivity corresponds to the self-inhibition  of the (four)
% neuronal populations constituting each region.
%
% Finally, we address the contribution of extrinsic and intrinsic
% pre-synaptic activity, laminar-specific contributions and the
% contributions of  inhibitory interneurons to the BOLD signal. This
% assessment uses Bayesian Model Reduction.
%
% 
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
 
% $Id: DEMO_dcm_fmri_nnm.m 6931 2016-11-16 12:09:58Z karl $
 
% tests of spatial models: 'ECD', 'LFP' or 'IMG'
%==========================================================================
try
    cd(fullfile(spm('Dir'),'tests','data','fMRI'))
catch
    try
        cd('C:\home\spm\DCM\DCM fMRI')
    catch
        cd('C:\Users\karl\Documents\SPM\DCM fMRI')
    end
end
load DCM_attention
 
DCM = rmfield(DCM,{'a','b','c','d','options'});
 
% spatial models
%==========================================================================
DCM.options.nmm = 'TFM';    % two regional populations (E and I)
 
 
% priors on connectivity
%--------------------------------------------------------------------------
DCM.a{1}        = [0 0 0;1 0 0;0 1 0];
DCM.a{2}        = [0 1 0;0 0 1;0 0 0];
DCM.b{1}(:,:,1) = [0 0 0;0 0 0;0 0 0];
DCM.b{1}(:,:,2) = [0 0 0;0 0 0;0 0 0];
DCM.b{1}(:,:,3) = [0 0 0;0 1 0;0 0 0];
DCM.b{2}(:,:,1) = [0 0 0;0 0 0;0 0 0];
DCM.b{2}(:,:,2) = [0 0 0;0 0 0;0 0 0];
DCM.b{2}(:,:,3) = [0 0 0;0 1 0;0 0 0];
DCM.c           = [1 0 0;0 1 0;0 0 1];
DCM.d           = [];
 
% Bayesian model inversion
%==========================================================================
DCM.options.maxit = 32;
DCM = spm_dcm_fmri_nmm(DCM);
 
% get posterior densities
%--------------------------------------------------------------------------
pE  = DCM.M.pE;
pC  = DCM.M.pC;
qE  = DCM.Ep;
qC  = DCM.Cp;
vC  = spm_unvec(diag(qC),qE);
 
 
% electrophysiological correlates
%==========================================================================
spm_figure('GetWin','Figure 1'); clf
 
% induced electrophysiological responses
%--------------------------------------------------------------------------
[y,lfp,csd,w] = spm_gen_fmri(DCM.Ep,DCM.M,DCM.U);
 
T   = (1:size(y,1))  *DCM.Y.dt;    % time of fMRI data
t   = (1:size(lfp,2))*DCM.U.dt;    % time of LFP data
j   = 4:64;                        % frequencies to report
n   = 2;                           % in the n-th region
ind = csd(:,:,n);                  % induced responses every TR
 
 
% posterior parameter estimates
%--------------------------------------------------------------------------
subplot(3,1,1), plot(T,y), hold on, plot(T,DCM.y + DCM.R,':'), hold off
title('Haemodynamic responses','fontsize',16), spm_axis tight
xlabel('Time (seconds)'),ylabel('Normalised BOLD')
 
subplot(3,1,2), plot(t,lfp')
title('Local field potential responses','fontsize',16), spm_axis tight
xlabel('Time (seconds)'),ylabel('depolarisation')
 
subplot(3,1,3), imagesc(T,w(j),spm_detrend(log(ind(:,j)))')
title('Induced  responses (log)','fontsize',16)
xlabel('Time (seconds)'),ylabel('Frequency (Hz)')
 
% principal spectral mode and haemodynamic response
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf
 
col = {'b','g','r'};
for i = 1:size(csd,3)
    ind     = spm_detrend(csd(:,:,i));
    [u,s,v] = spm_svd(ind);
    f       = sign(u'*spm_detrend(y(:,1)));
    u       = u*diag(f);
    v       = v*diag(f);
    str     = sprintf('First eigenmode: %s',DCM.xY(i).name);
    subplot(3,2,2*i - 1), plot(w,v(:,1)), axis square, spm_axis tight, box off
    title(str,'fontsize',16), xlabel('Frequency'),ylabel('Weight')
    subplot(3,2,2*i - 0), plot(u(:,1),DCM.Y.y(:,i),['.' col{i}],'Markersize',8), hold on
    subplot(3,2,2*i - 0), plot(u(:,1),DCM.Y.y(:,i),[':' col{i}]), hold off
    title('Haemodynamic correlates','fontsize',16), axis square
    xlabel('Eigenvariate'),ylabel('Normalised BOLD')
end
 
% Bayesian model comparison and reduction
%==========================================================================
 
% neurovascular coupling parameters
%--------------------------------------------------------------------------
j        = spm_fieldindices(qE,'J');
dcm.M.pE = DCM.M.pE.J;
dcm.M.pC = DCM.M.pC.J;
dcm.Ep   = DCM.Ep.J;
dcm.Cp   = DCM.Cp(j,j);
 
% posterior parameter estimates
%--------------------------------------------------------------------------
spm_dcm_bmr_all(dcm)
 
subplot(3,2,5), imagesc(spm_cov2corr(qC(j,j))), axis square
title('Neurovascular coupling','fontsize',16), 
xlabel('pre-synaptic input'),ylabel('Contribution')
 
 
% The effects of attention
%==========================================================================
spm_figure('GetWin','Figure 3'); clf
 
% extrinsic verses intrinsic and superficial versus deep
%--------------------------------------------------------------------------
i  = spm_fieldindices(qE,'N');
j  = spm_find_pC(vC.N,qE.N,'B');
j  = i(j);
 
% posterior parameter estimates
%--------------------------------------------------------------------------
subplot(2,2,1), spm_plot_ci(qE,qC,[],j), axis square
title('Attentional modulation','fontsize',16), 
xlabel('Connection'),ylabel('Effect size (log scaling)')
 
% Bayesian model comparison
%--------------------------------------------------------------------------
clear rE rC F
 
rE    = pE.N.B; 
rC{1} = pC.N.B;                                           % full model
rC{2} = pC.N.B; rC{2}{1}(:,:,3) = 0;                      % superficial
rC{3} = pC.N.B; rC{3}{2}(:,:,3) = 0;                      % deep
rC{4} = pC.N.B; rC{4}{1}(:,:,3) = 0; rC{4}{2}(:,:,3) = 0; % both
 
for i = 1:numel(rC)
    F(i,1) = spm_log_evidence(qE.N.B,vC.N.B,rE,rC{1},rE,rC{i});
end
p = spm_softmax(F);
 
% show results
%--------------------------------------------------------------------------
subplot(2,2,2), bar(p), axis square
title('Model comparison','fontsize',16), 
xlabel('Model'),ylabel('Posterior probability')
set(gca,'XTickLabel',{'full','supra','infra','null'})
