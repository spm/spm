% Demo for a bird songs: In this example, we show that DEM can not only
% estimate the hidden states of an autonomous system but can also
% deconvolve dynamics changes in its control parameters.  We illustrate
% this using a slow Lorentz attractor to drive a fast one; both  showing
% deterministic chaos.  We endow the simulations with a little ethological
% validity by using the states of the fast Lorentz attractor as control
% variables in the syrinx of a song bird (usually these would control a van
% der Pol oscillator model). We will look at the true and inferred songs
% with and without the last part missing.  When sonograms are displayed the
% song can be played by a mouse click on the image. The final plots show
% simulated event related potential to show that there is a marked
% responses (prediction error) of the system when an expected ‘syllable’ is
% omitted. This demonstrates the implicit sequence-decoding of input
% streams, using uncontrollable state-space models
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_song_inference.m 3655 2009-12-23 20:15:34Z karl $
 
 
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
% P(1): Prandtl number
% P(2): 8/3
% P(3): Rayleigh number
 
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
