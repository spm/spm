function LAP = DEM_demo_texture
% This demonstration considers the figure-ground segregation problem where, 
% crucially, a figure is defined texturally – in terms of its second order 
% statistics; in other words, a visual object is manifest in terms of its 
% texture or spectral power density in the spatial domain. This definition 
% precludes any first-order attributes; such as increased luminance. This 
% sort of problem is common in the inverse literature and is usually solved 
% using a prior on the [co]variance of random fluctuations generating data. 
% Here, we simulate a contiguous object, whose texture is determined by the 
% variance of random fluctuations in luminance – and the variance (or 
% precision) is modulated by Gaussian basis functions. The resulting signal 
% is mixed with uniform Gaussian noise to produce sensory data. These 
% (one-dimensional) data are then subject to Bayesian inversion using 
% generalized predictive coding – (as implemented in spm_LAP) – to recover 
% the underlying object.
% 
% Technically, this scheme optimizes expectations of the hidden causes of 
% the data, which are the amplitude of radial basis functions controlling 
% the precision of retinotopic signals. Heuristically, the figure is 
% recognized by selectively attending to sensory input from the figure. 
% This enables sensory noise to be suppressed in unattended parts of the 
% visual field. However, this form of attention is distinct from simply 
% boosting sensory precision (the precision of sensory prediction errors) 
% as in simulations of the Posner paradigm or biased competition. Here, 
% the hidden causes are optimized in a way that renders them less precise 
% and therefore more sensitive to ascending (prediction error) sensory 
% input. This illustrates the functional importance of the relative 
% precision of sensory and extrasensory prediction errors in modulating 
% the influence of ascending sensory information that competes to influence 
% posterior expectations.
%
% PS: for a 2-D simulation delete 'return' below.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_texture.m 6235 2014-10-12 10:03:05Z karl $
 

% Create a generative model:
%==========================================================================
rng('default')
                                       
% level 1; textured stimulus with noise (log precision of four)
%--------------------------------------------------------------------------
G(1).v  = zeros(128,1);                   % output channels (stimuli)
G(1).V  = 16;                             % error precision (noise)
G(1).g  = inline('spm_conv(v,2)','x','v','P');


% level 2; underlying causes (three Gaussian patches)
%--------------------------------------------------------------------------
G(2).v  = zeros(128,1);                   % textured stimulus
G(2).V  = [];
G(2).ph = @ph1;
G(2).g  = inline('zeros(128,1)','x','v','P');


% level 2; amplitude and size
%--------------------------------------------------------------------------
G(3).v  = zeros(3,1);
G(3).V  = exp(8);

% evaluate G to generate stimulus
%--------------------------------------------------------------------------
U       = [1;1;0]*(8 + 4);                % amplitude and size
LAP     = spm_DEM_generate(G,U);


 
% invert to simulate  predictive coding
%==========================================================================
LAP.M(3).v = U - U + 8;                           % use correct starting estimates
LAP.M(3).V = 1/2;
LAP        = spm_LAP(LAP);

% plot results
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');
spm_DEM_qU(LAP.qU,LAP.pU)

subplot(3,2,1), title('prediction and error','FontSize',16)
subplot(3,2,2), title('signal plus noise','FontSize',16)
subplot(3,2,3), title('posterior beliefs','FontSize',16)
subplot(3,2,4), title('hidden cause (1st order)','FontSize',16)
subplot(3,2,5), title('posterior beliefs','FontSize',16)
subplot(3,2,6), title('hidden cause (2nd order)','FontSize',16)

f  = exp(-ph1([],LAP.qU.v{3},[],LAP.M)/2);
subplot(3,2,3), hold on, plot(f,'r')
f  = LAP.G(1).g([],LAP.pU.v{2},[]);
subplot(3,2,1), hold on, plot(f,'b:')
f  = exp(-ph1([],LAP.pU.v{3},[],LAP.G)/2);
subplot(3,2,4), hold on, plot(f,'r')


return
 
% Create a generative model: two dimensions
%==========================================================================
clear 
rng('default')
                                       
% level 1; textured stimulus with noise (log precision of four)
%--------------------------------------------------------------------------
G(1).v  = zeros(32^2,1);                   % output channels (stimuli)
G(1).V  = 16;                              % error variances (noise)
G(1).g  = inline('spm_conv(v,2)','x','v','P');


% level 2; underlying causes (three Gaussian patches)
%--------------------------------------------------------------------------
G(2).v  = zeros(32^2,1);                   % textured stimulus
G(2).V  = [];
G(2).ph = @ph2;
G(2).g  = inline('zeros(32^2,1)','x','v','P');


% level 2; amplitude and size
%--------------------------------------------------------------------------
G(3).v  = zeros(9,1);
G(3).V  = exp(8);


% evaluate G to generate stimulus
%--------------------------------------------------------------------------
U       = spm_vec([0 0 0; 1 0 0; 1 1 0]*8);
LAP     = spm_DEM_generate(G,U);
 
% invert to simulate  predictive coding
%==========================================================================
LAP.M(3).v = U;
LAP.M(3).V = 1/8;
LAP        = spm_LAP(LAP);


% plot results
%--------------------------------------------------------------------------
spm_DEM_qU(LAP.qU,LAP.pU)
spm_figure('GetWin','Figure 2');
spm_DEM_qU(LAP.qU,LAP.pU)

qp  = exp(-ph2([],LAP.qU.v{3},[],LAP.M)/2);
pp  = exp(-ph2([],LAP.pU.v{3},[],LAP.G)/2);
qf  = LAP.M(1).g([],LAP.qU.v{2},[]);
pf  = LAP.G(1).g([],LAP.pU.v{2},[]);
n   = 32;

subplot(3,3,1), imagesc(reshape(pf,32,32)); axis square, title('feature','FontSize',16)
hold on, plot([0 n],([n n] + 1)/2,':w',([n n] + 1)/2,[0 n],':w'), hold off
subplot(3,3,2), imagesc(reshape(LAP.Y,32,32));       axis square, title('stimulus','FontSize',16)
hold on, plot([0 n],([n n] + 1)/2,':w',([n n] + 1)/2,[0 n],':w'), hold off
subplot(3,3,3), imagesc(reshape(qf,32,32)); axis square, title('percept','FontSize',16)
hold on, plot([0 n],([n n] + 1)/2,':w',([n n] + 1)/2,[0 n],':w'), hold off
subplot(3,2,3), imagesc(reshape(qp,32,32)); axis square, title('predicted variance','FontSize',16)
hold on, plot([0 n],([n n] + 1)/2,':w',([n n] + 1)/2,[0 n],':w'), hold off
subplot(3,2,4), imagesc(reshape(pp,32,32)); axis square, title('true variance','FontSize',16)
hold on, plot([0 n],([n n] + 1)/2,':w',([n n] + 1)/2,[0 n],':w'), hold off
subplot(3,2,6), title('true causes','FontSize',16)

return


function p = ph1(x,v,h,M)
% returns log precision
%__________________________________________________________________________
n     = numel(v);
x     = 1:128;
c     = linspace(48,80,n);
p     = x - x;
for i = 1:n
    b = exp(-((x - c(i)).^2)/128);
    p = p + (b - mean(b))*v(i);
end

p  = 8 - p(:);

return


function p = ph2(x,v,h,M)
% returns log precision
%__________________________________________________________________________
n       = sqrt(numel(v));
v       = reshape(v,n,n);
x       = 1:32;
[x1,x2] = ndgrid(x);
c       = linspace(8,24,n);
p       = x1 - x1;

for i = 1:size(v,1)
    for j = 1:size(v,1)
        b = exp(-((x1 - c(i)).^2 + (x2 - c(j)).^2)/32);
        p = p + (b - mean(mean(b)))*v(i,j);
    end
end

p  = 8 - p(:);





