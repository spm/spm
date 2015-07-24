function DEM_demo_ontology
% Demo for a bird songs: In this example, we simulate local field potential
% using the prediction error from the song-bird example below. We look at
% these responses under natural stimuli and after removing the second
% level of the hierarchy to show it is necessary for veridical perception.
% We then repeat but omitting dynamical priors by forsaking generalised
% coordinates
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_demo_ontology.m 6506 2015-07-24 10:26:51Z karl $


% hierarchical non-linear generative model
%==========================================================================
rng('default')

% length of trajectory
%--------------------------------------------------------------------------
N        = 64;                      % length of stimulus (bins)

% serial correlations among latent sstates
%--------------------------------------------------------------------------
M(1).E.s  = 1/2;
M(1).E.n  = 2;
M(1).E.K  = 1/16;
M(1).E.nE = 32;

% therapeutic intervention
%--------------------------------------------------------------------------
U        = spm_phi(((1:N) - N/2));

% level 1
%--------------------------------------------------------------------------
P.A      = randn(4,2)/32;
P.B      = [0 0 1 1; 0 1 0 1];
M(1).g   = @(x,v,P)[spm_phi(P.A*exp(v)); spm_softmax(-sum((P.B - v*ones(1,4)).^2)')];
M(1).pE  = P;
M(1).V   = exp(8);


% level 2
%--------------------------------------------------------------------------
P.A      = [10 24 1];
P.B      = [2 0;1/2 1];
x        = [2; 4; 32];
M(2).f   = @(x,v,P)[-P.A(1) P.A(1) 0; ((1 - v*P.A(3))*P.A(2) - x(3)) -1 0; x(2) 0 -8/3]*x/32;
M(2).g   = @(x,v,P) P.B*x([2 3])/16;
M(2).x   = x;
M(2).pE  = P;
M(2).pC  = 0;
M(2).V   = exp(8);
M(2).W   = exp(4);

% level 3
%--------------------------------------------------------------------------
M(3).v   = U(:,1);
M(3).V   = exp(16);


% illustrate trajectories with and without therapy
%==========================================================================

% set subject specific parameters
%--------------------------------------------------------------------------
P      = {M.pE};
P{2}.A = P{2}.A + [0 8 0];


% natural progression without therapy
%--------------------------------------------------------------------------
DEM.U  = U*0;
DEM    = spm_DEM_generate(M,DEM.U,P);
DEM    = spm_DEM(DEM);

spm_figure('GetWin','Figure 1'); clf
spm_DEM_plot(DEM)

% repeat with therapy
%--------------------------------------------------------------------------
DEM.U  = U;
DEM    = spm_DEM_generate(M,DEM.U,P);
DEM    = spm_DEM(DEM);

spm_figure('GetWin','Figure 2'); clf
spm_DEM_plot(DEM)


% prognosis and prediction
%==========================================================================

% Estimate subject specific parameters at presentation
%--------------------------------------------------------------------------
PEM         = DEM;
PEM.U       = PEM.U(:,1:N/2);
PEM.Y       = PEM.Y(:,1:N/2);
PEM.M(2).pC = diag([1 1 1 0 0 0 0]);
PEM         = spm_DEM(PEM);

spm_figure('GetWin','Figure 3'); clf
spm_DEM_qP(PEM.qP,PEM.pP)

% update priors and predict with zero precision data
%--------------------------------------------------------------------------
PEM         = spm_ADEM_update(PEM,0);
PEM.M(1).V  = 0;
PEM.M(2).W  = exp(8);
PEM.M(2).pC = 0;

PEM.U       = U(N/4:end);
PEM.Y       = zeros(8,size(PEM.U,2));
PEM         = spm_DEM(PEM);

spm_figure('GetWin','Figure 4'); clf
spm_DEM_plot(PEM)

% replear with out therapeutic intervention
%--------------------------------------------------------------------------
PEM.U = U(N/4:end)*0;
PEM   = spm_DEM(PEM);

spm_figure('GetWin','Figure 5'); clf
spm_DEM_plot(PEM)

u     = linspace(0,2,8);
for i = 1:length(u)
    
    % ppredict a response to treatment level
    %----------------------------------------------------------------------
    PEM.U   = U(N/4:end)*u(i);
    PEM     = spm_DEM(PEM);
    
    % plot response trajectory
    %----------------------------------------------------------------------
    spm_figure('GetWin','Figure 5'); subplot(3,2,6), hold on
    plot(PEM.qU.v{2}(1,:),PEM.qU.v{2}(2,:),'r:')
    
    % record predictive efficacy
    %----------------------------------------------------------------------
    R(i)    = PEM.qU.v{1}(5,end);
    
end

spm_figure('GetWin','Figure 5'); subplot(3,2,5), hold off
bar(u,R)
title('dose response curve','fontsize',16)
xlabel('level of therapeutic intervention')
xlabel('probability of remission')

% return





% model identification and selection
%==========================================================================
n       = 8;                 % nnumber of subjects
M(2).pC = 1;
for i = 1:n
    
    % set subject specific parameters
    %----------------------------------------------------------------------
    P      = {M.pE};
    P{2}.A = P{2}.A + [0 8 0] + randn(1,3);

    % repeat with therapy
    %----------------------------------------------------------------------
    DCM{i,1}   = spm_DEM_generate(M,U,P);
    DCM{i,1}.U = U;
    
end

for i = 1:4
    DCM       = spm_dcm_fit(DCM);
    [PEB,DCM] = spm_dcm_peb(DCM);
end

for i = 1:n
    
    % set subject specific parameters
    %----------------------------------------------------------------------
    qA(:,i) = DCM{i}.qP.P{2}.A(:);
    pA(:,i) = DCM{i}.pP.P{2}.A(:);
    qB(:,i) = DCM{i}.qP.P{2}.B(:);
    pB(:,i) = DCM{i}.pP.P{2}.B(:);

end

PEB.Ep  = spm_unvec(PEB.Ep,DCM{1}.M(2).pE)
spm_dcm_bmr_all(PEB,{'B'})


% show songs and prediction error (ERP)
%==========================================================================
function spm_DEM_plot(DEM)

% plot hidden states and causes
%--------------------------------------------------------------------------
if isfield(DEM,'pU')
    spm_DEM_qU(DEM.qU,DEM.pU)
else
    spm_DEM_qU(DEM.qU)
end

% supplement with trajectories
%--------------------------------------------------------------------------
subplot(6,2,2), imagesc(DEM.qU.v{1}(1:4,:))
title('symptom profile and diagnosis','fontsize',16)
subplot(6,2,4), imagesc(DEM.qU.v{1}(5:8,:))

% supplement with trajectories
%--------------------------------------------------------------------------
a   = [-1 3 -1 3];
vi  = linspace(a(1),a(2),64);
vj  = linspace(a(3),a(4),64);
for i = 1:length(vi)
    for j = 1:length(vj)
        x      = DEM.M(1).x;
        v      = [vi(i); vj(j)];
        p      = DEM.M(1).g(x,v,DEM.qP.P{1});
        p      = p(5:8);
        s(j,i) = p'*log(p);
        [p,q]  = max(p);
        d(j,i) = q;
    end
end
d    = d.*(min(s(:)) - s);

subplot(3,2,6), imagesc(a(1:2),a(3:4),d), axis xy, hold on
plot(DEM.qU.v{2}(1,:),DEM.qU.v{2}(2,:),'r'), hold off
title('latent psychopathology','fontsize',16)


subplot(3,2,4),title('latent pathophysiology','fontsize',16)
subplot(3,2,3),title('latent psychopathology','fontsize',16)
subplot(3,2,1),title('predicted symptoms and error','fontsize',16)
subplot(3,2,5),title('pathology and therapy','fontsize',16)
set(gca,'YLim',[-.2 1.2])

