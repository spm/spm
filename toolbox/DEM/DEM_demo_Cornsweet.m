function DEM_demo_Cornsweet
% The Cornsweet effect: This demo illustrates the inference underlying the
% Cornsweet effect or illusion. It exploits formal priors on the spatial
% contiguity of the illuminant and reflectance; where the illuminant does not
% have edges, but the reflectance can. This is implemented using a
% discrete cosine set (DCT) as the spatial basis for the illuminant and a 
% (Haar) Discrete Wavelet transform (DWT) for the reflectance. Appropriate
% shrinkage priors on the (implicit) transform coefficients ensure that the
% explanation for visual input (reflectance times illuminant) assigns edges
% to the reflectance; thereby producing the Cornsweet effect.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_Cornsweet.m 4851 2012-08-20 15:03:48Z karl $
 
 
% Illustrate the Cornsweet effect
%==========================================================================
spm_figure('GetWin','Figure 1');
colormap((1:255)'*[1 1 1]/255)
 
% basic profile
%--------------------------------------------------------------------------
nx    = 64;
YX    = [zeros(1,nx*3/8) 1:(nx/8) -(nx/8):-1 zeros(1,nx*3/8)]*8/nx;
YT    = ones(32,1);
 
% present with different contrasts (C)
%--------------------------------------------------------------------------
C     = [1 4 32 64];
nc    = length(C);
for i = 1:length(C)
 
    subplot(nc,2,2*i - 1)
    image(YT*YX*C(i) + 128)
    title('Cornsweet effect','FontSize',16)
    xlabel('eccentricity ','FontSize',12)
    axis square
 
    subplot(nc,2,2*i)
    plot(YX*C(i))
    title('Reflectance ','FontSize',16)
    xlabel('eccentricity ','FontSize',12)
    axis square, axis([1 nx -64 64]);
    
end
 
 
% Simulations
%==========================================================================
 
 
% Basis functions of illuminant and reflectance
%--------------------------------------------------------------------------
nx      = 32;                               % number of pixels
P.R     = spm_dwtmtx(nx,1,1);               % DWT for hidden causes
P.I     = spm_dctmtx(nx,3);                 % DCT for hidden causes
P.I     = P.I/diag(max(P.I));
nr      = size(P.R,2);
ni      = size(P.I,2);
W       = log2(nx./sum(P.R > 0))*4;         % scale of reflectance
W(1)    = 16;
 
 
% initial hidden states and causes
%--------------------------------------------------------------------------
x       = sparse(nr,1);
v       = sparse(ni + nr,1);
v(1)    = -1;
ii      = (1:ni);
ir      = (1:nr) + ni;
 
% level 1
%--------------------------------------------------------------------------
M(1).f  = inline('(v(4:end) - x)/16','x','v','P');
M(1).g  = inline('exp(P.R*x + P.I*v(1:3))','x','v','P');
M(1).pE = P;                                % The prior expectation
M(1).x  = x;                                % The prior expectation
M(1).V  = exp(6);                           % error precision (data)
M(1).W  = exp(10);                          % error precision (motion)
M(1).xP = diag(exp(W));                     % error precision (motion)
 
% level 2
%--------------------------------------------------------------------------
M(2).v  = v;                                % hidden causes 
M(2).V  = exp(0);                           % error precision (cause)
 
% Create stimulus
%==========================================================================
 
% Create stimulus (with spatial and temporal envelopes YX and YT)
%--------------------------------------------------------------------------
M(1).E.n = 1;
N        = 64;                             % length of sequence
t        = ([1:N] - N/4)*8;                % time (ms)
YX       = [zeros(1,nx*3/8) 1:(nx/8) -(nx/8):-1 zeros(1,nx*3/8)]*8/nx;
YX       = YX'/4 + 1;
YX       = exp(P.R*pinv(P.R)*log(YX));
YT       = exp((tanh(([1:N] - N/3)/8) - 1)*4);
YT       = exp(-([1:N] - N/2).^2/(2*(N/8)^2));
Y        = YX*YT;
 
% invert
%==========================================================================
DEM.M    = M;
DEM.Y    = Y;
DEM      = spm_DEM(DEM);
 
 
% render true and perceived stimuli
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf
 
subplot(2,2,1)
PI  = P.I*DEM.qU.v{2}(ii,:);
imagesc(exp(PI))
title('Perceived illunimant','FontSize',16)
xlabel('time (bins)','FontSize',12)
ylabel('eccentricity ','FontSize',12)
axis square
 
subplot(2,2,2)
PR = P.R*DEM.qU.x{1};
imagesc(exp(PR))
title('Perceived reflectance','FontSize',16)
xlabel('time (bins)','FontSize',12)
ylabel('eccentricity ','FontSize',12)
axis square
 
subplot(2,2,3)
imagesc(DEM.qU.v{1})
title('Predicted stimulus','FontSize',16)
xlabel('time (bins)','FontSize',12)
ylabel('eccentricity ','FontSize',12)
axis square
 
 
% first order effect: (R1) - Cornsweet
%--------------------------------------------------------------------------
C     = sparse(1,nx);
i     = nx/2 - nx/4 - 1;
C(i)  =  1;
i     = nx/2 + nx/4 + 2;
C(i)  = -1;
R1    = C*P.R;
 
% second order effect: (R2) - Mach bands
%--------------------------------------------------------------------------
i     = nx/2 - nx/8 - 1;
C(i)  = -1;
i     = nx/2 + nx/8 + 2;
C(i)  =  1;
R2    = C*P.R;
 
D1    = R1*DEM.qU.x{1};
D2    = R2*DEM.qU.x{1};
for i = 1:N
    V1(i) = R1*DEM.qU.S{i}*R1';
    V2(i) = R2*DEM.qU.S{i}*R2';
end
 
 
subplot(2,2,4)
spm_plot_ci(D1,V1,t), hold on
spm_plot_ci(D2,V2,t), hold off
title('Perceived reflectance','FontSize',16)
xlabel('time (ms)','FontSize',12)
ylabel('eccentricity ','FontSize',12)
axis square
spm_axis tight
drawnow
 
 
% Cycle over different levels of visual precision (cf contrast)
%==========================================================================
C   = -2:2:12;
for i = 1:length(C)
    
    % Change precision (contrast) and invert
    %----------------------------------------------------------------------
    SIM{i}        = DEM;
    SIM{i}.M(1).V = exp(C(i));
    SIM{i}        = spm_DEM(SIM{i});
 
    % Record conditional peripheral reflectance difference
    %----------------------------------------------------------------------
    d1(i) = R1*SIM{i}.qU.x{1}(:,N/2);
    v1(i) = R1*SIM{i}.qU.S{N/2}*R1';
    d2(i) = R2*SIM{i}.qU.x{1}(:,N/2);
    v2(i) = R2*SIM{i}.qU.S{N/2}*R2';
 
end
 

 
% show Cornsweet effect as a function of precision (contrast)
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3'); clf
 
subplot(3,1,1)
spm_plot_ci(d1,v1,C),                hold on
plot(C,d1,'ob',C,C*0,'LineWidth',2), hold off
title('Conditional difference (Cornsweet)','FontSize',16)
xlabel('log-precision (contrast)','FontSize',12)
ylabel('reflectance difference','FontSize',12)
spm_axis tight square
 
subplot(3,1,2)
spm_plot_ci(d2,v2,C),                hold on
plot(C,d2,'or',C,C*0,'LineWidth',2), hold off
title('Conditional difference (Mach Band)','FontSize',16)
xlabel('log-precision (contrast)','FontSize',12)
ylabel('reflectance difference','FontSize',12)
spm_axis tight square
 
% and associated percepts
%--------------------------------------------------------------------------
colormap([1:255]'*[1 1 1]/255)
j     = [1 5 7];
nj    = length(j);
for i = 1:nj
   
    % predictions
    %----------------------------------------------------------------------
    PR    = P.R*SIM{j(i)}.qU.x{1};
    PR    = PR*256 + 128;
    
    subplot(3,nj,i + 2*nj)
    image(PR)
    ylabel('eccentricity','FontSize',12)
    axis square
 
end
 
% and associated prediction errors (cf ERPs)
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 4'); clf
as    = 0;
ax    = 0;
av    = 0;
for i = 1:nj
    
    
    % predictions and errors
    %----------------------------------------------------------------------
    Ev = SIM{j(i)}.qU.z{2}(ir,:);
    Es = SIM{j(i)}.qU.z{1};
    av = max(max(abs(Ev(:))),av);
    as = max(max(abs(Es(:))),as);
    
    subplot(3,nj,i)
    plot(t,Es,'r')
    title('Error (sensory)','FontSize',14)
    xlabel('time (ms)','FontSize',12)
    ylabel('error (state)','FontSize',12)
 
    subplot(3,nj,i + nj)
    plot(t,Ev,'r')
    title('Error (hidden states)','FontSize',14)
    xlabel('time (ms)','FontSize',12)
    ylabel('error (state)','FontSize',12)
    
 
end
 
 
% adjust axes
%--------------------------------------------------------------------------
for i = 1:nj
    
    subplot(3,nj,i)
    axis square, axis([t(1) t(end) -as*1.2 as*1.2])
    
    subplot(3,nj,i + nj)
    axis square, axis([t(1) t(end) -av*1.2 av*1.2])
 
end
 
 
% redraw inference for effect
%--------------------------------------------------------------------------
spm_figure('GetWin','DEM');
spm_DEM_qU(SIM{7}.qU)

 
% Fitting behavioural data
%==========================================================================
try, load cornsweet_data, catch, return, end

% save simulation results
%--------------------------------------------------------------------------
sim.Con    = C;               % simulated contrast
sim.Ecs    = d1;              % conditional expectation of Cornsweet effect
sim.Emb    = d2;              % conditional expectation of mach band effect
sim.Vcs    = v1;              % conditional dispersion of Cornsweet effect
sim.Vmb    = v2;              % conditional dispersion of mach band effect
 
% add emprical contrasts that were used
%--------------------------------------------------------------------------
sim.Con_cs = cornsweet.cornsweet;
sim.Con_mb = mach.machContrast;
G.sim      = sim;     % place in model
 
% empirical responses
%--------------------------------------------------------------------------
Y          = struct;
Y.y{1}     = cornsweet.stepMatch(:);
Y.y{2}     = mach.pSeeMach(:);
Y.Q        = spm_Ce([length(Y.y{1}) length(Y.y{2})]); % error precisions
 
 
% model parameters
%--------------------------------------------------------------------------
B.thresh = -2;                % log-contrast threshold for Mach band
B.sig    = 1;                 % parameter of contrast function
B.off    = 2;                 % parameter of contrast function
B.con    = 1/32;              % scaling of reported contrast
B.sen    = 4;                 % parameter of response probability
nb       = length(spm_vec(B));
 
 
% setup model of empirical responses
%--------------------------------------------------------------------------
G.IS = 'spm_cornsweet';       % function name f(P,M,U) - generative model
G.pE = B;                     % prior expectation of model parameters
G.pC = speye(nb,nb)/8;        % prior covariance  of model parameters
G.hE = [16; 8];               % prior expectation of log-precisions
G.hC = exp(-4);               % prior covariance of log-precisions
 
 
% invert model of empirical responses and plot
%--------------------------------------------------------------------------
Ep   = spm_nlsi_GN(G,[],Y);
spm_cornsweet(Ep,G,Y);




return

% plot behavioural results alone.
%--------------------------------------------------------------------------
subplot(2,2,1)
SE   = cornsweet.error;
errorbar(cornsweet.cornsweet,cornsweet.stepMatch,SE/2)
xlabel('empirical contrast','Fontsize',12)
ylabel('reported contrast','Fontsize',12)
title('Cornsweet','Fontsize',16)
set(gca,'XLim',[0 .15]);
axis square

subplot(2,2,2)
SE   = mach.error;
errorbar(mach.machContrast,mach.pSeeMach,SE/2)
xlabel('empirical contrast','Fontsize',12)
ylabel('report probability','Fontsize',12)
title('Mach Bands','Fontsize',16)
set(gca,'XLim',[0 .15]);
axis square
 
% Ancillary code for producing figures:
%==========================================================================
spm_figure('GetWin','paper'); clf
 
n = 128;
a = 4;
i = tanh(([1:n] - n/2)/(n/8))/a;
r = kron([1 -1],ones(1,n/2))/a;
s = exp(i + r);
c = zeros(1,n);
 
subplot(5,1,1), plot(i),     spm_axis tight, title('illuminant')
subplot(5,1,2), plot(r),     spm_axis tight, title('reflectant')
subplot(5,1,3), plot(s),     spm_axis tight, title('stimulus')
subplot(5,1,4), plot(i + r), spm_axis tight, title('reflectant')
subplot(5,1,5), plot(c),     spm_axis tight, title('illuminant')
 
 
% Illustrate the generative model (for 1D)
%==========================================================================
spm_figure('GetWin','paper'); clf
 
DEM      = spm_DEM_generate(M,64,P,{16 0},{16});
 
subplot(3,2,1)
PI  = P.I*DEM.pU.v{2}(ii,:);
imagesc(exp(PI))
title('Illunimant','FontSize',16)
xlabel('time (bins)','FontSize',12)
ylabel('location ','FontSize',12)
axis square
 
subplot(3,2,3)
PR = P.R*DEM.pU.x{1};
imagesc(exp(PR))
title('Reflectant','FontSize',16)
xlabel('time (bins)','FontSize',12)
ylabel('location ','FontSize',12)
axis square
 
subplot(3,2,5)
imagesc(DEM.pU.v{1})
title('Stimulus','FontSize',16)
xlabel('time (bins)','FontSize',12)
ylabel('location ','FontSize',12)
axis square
 
 
subplot(3,2,2), imagesc(P.I), axis image off
subplot(3,2,4), imagesc(P.R), axis image off
 
 
% Illustrate the generative model (for 2D)
%==========================================================================
nx    = 128;
P.R   = spm_dwtmtx(nx,2);                 % DWT for hidden causes
P.I   = spm_dctmtx(nx,3);                 % DCT for hidden causes
I     = 0;
R     = 0;
for i = 1:size(P.I,2)
    for j = 1:size(P.I,2)
        I = I + (P.I(:,i)*P.I(:,j)')*randn;
    end
end
for i = 1:size(P.R,2)
    for j = 1:size(P.R,2)
        v = sum(P.R(:,i) > 0)*sum(P.R(:,j) > 0);
        R = R + (P.R(:,i)*P.R(:,j)')*randn*sqrt(v);
    end
end
 
R   = R/std(R(:))/2;
I   = I/std(I(:))/2;
 
% Plot
%--------------------------------------------------------------------------
spm_figure('GetWin','Paper');
 
subplot(2,2,1)
imagesc(exp(I))
title('Illuminant ','FontSize',16)
axis square
 
subplot(2,2,2)
imagesc(exp(R))
title('Reflectance ','FontSize',16)
axis square
 
subplot(2,1,2)
imagesc(exp(I + R))
title('Image','FontSize',16)
axis square, drawnow


