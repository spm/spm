% Perceptual categorisation of bird songs: The generative model of 
% birdsong used in this simulation comprises a Lorenz attractor with two 
% control parameters (or hidden causes), which, in turn, delivers two 
% control parameters to a synthetic syrinx to produce ‘chirps’ that are 
% modulated in amplitude and frequency.  The chirps were then presented 
% as a stimulus to a synthetic bird to see if it could infer the 
% underlying causal states and thereby categorise the song. This entails 
% minimising free energy by changing the internal representation of the 
% control parameters. Each simulated song comprises a series of chirps 
% whose frequency and number fall progressively from song a to song c, 
% as a causal state (known as the Raleigh number) is decreased.  The 
% simulations show that the causes are identified after about 600 
% milliseconds with high conditional precision (90% confidence intervals 
% are shown in grey). These simulations illustrate the nature of 
% perceptual categorisation under generalised predictive coding: Here, 
% recognition corresponds to mapping from a continuously changing and 
% chaotic sensory input to a fixed point in perceptual space.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_song_inference.m 4170 2011-01-24 18:37:42Z karl $
 
 
% Hierarchical non-linear generative model (dynamic & chaotic)
%==========================================================================
spm_figure('Getwin','Graphics');
clear M
 
% timing
%--------------------------------------------------------------------------
N        = 64;                       % length of stimulus (bins)
dt       = 1/64;                     % time bin (seconds)
t        = [1:N]*dt;

% correlations
%--------------------------------------------------------------------------
M(1).E.s = 1;
M(1).E.K = exp(-2);
 
% level 1
%--------------------------------------------------------------------------
% Prandtl number  = 10
% P(2)            = v(2)
% Rayleigh number = v(1) - 4: 
 
x       = [0.9; 0.8; 2];
M(1).f  = ' [-10 10 0; (v(1) - 4 - x(3)) -1 0; x(2) 0 -v(2)]*x/16;';
M(1).g  = 'x([2 3])';
M(1).x  = x;
M(1).V  = exp(4);
M(1).W  = exp(4);
 
 

% level 3
%--------------------------------------------------------------------------
M(2).v  = [0 0]';
M(2).V  = exp(-4);
 

% create data and invert three songs
%==========================================================================
S     = spm_phi(([1:N] - N/8)/(N/32));
P     = [32 26 16;
         2 8/3 8/3];
str = {'Song A','Song B','Song C'};

for i = 1:size(P,2)

    % create innovations & add causes
    %----------------------------------------------------------------------
    U(1,:)   = P(1,i)*S;
    U(2,:)   = P(2,i)*S;
    DEM      = spm_DEM_generate(M,U,{[] [] []},{4 16 16},{16 16 []});


    % DEM estimation and display
    %======================================================================
    DEM      = spm_DEM(DEM);

    % show song
    %----------------------------------------------------------------------
    spm_DEM_qU(DEM.qU,DEM.pU)
    subplot(2,2,4),colormap('pink')
    spm_DEM_play_song(DEM.qU,N*dt);
    axis square
    title('percept','Fontsize',16)  
     
    % record song
    %----------------------------------------------------------------------
    spm_figure('Getwin','Graphics');
    if i == 1; clf; end
    
    subplot(3,size(P,2),i)
    spm_DEM_play_song(DEM.qU,N*dt);
    axis square
    title(str{i},'Fontsize',16)  
    
    subplot(3,1,2)
    spm_plot_ci(t,DEM.qU.v{2},DEM.qU.C)
    text(t(end) + 1/8,DEM.qU.v{2}(1,end - 8),str{i},'Fontsize',16)

    axis square
    hold on
    
    subplot(3,1,3)
    spm_plot_ci(t,DEM.qU.v{2}(:,end - 8),DEM.qU.C(end - 8))
    axis square
    hold on
    plot(P(1,i),P(2,i),'or')
    text(P(1,i),P(2,i) + 1/2,str{i},'Fontsize',16)

end

drawnow, disp(' '),disp('Click sonograms to play songs'),disp(' ')
