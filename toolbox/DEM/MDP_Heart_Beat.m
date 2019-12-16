function MDP = MDP_Heart_Beat
% This simulation uses a Markov Decision process formulation of active
% inference to demonstrate the interaction between interoceptive and
% exteroceptive perception. This relies upon the fact that the function of
% exteroceptive sense organs depends upon oscillatory cycles in
% interoceptive states. The example used here is the change in retinal
% blood flow, and its influence on vision, during a cardiac cycle.

rng default

% Precisions
%--------------------------------------------------------------------------
wI = 1; % Interoceptive (inverse) volatility
wE = 1; % Exteroceptive (inverse) volatility
xi = 3; % Preferences
zI = 0.9; % Interoceptive sensory precision
zE = 0.9; % Exteroceptive sensory precision

% Prior beliefs about initial states
%--------------------------------------------------------------------------
D{1} = [1 1 1]';     % Interoceptive states (diastole, diastole, systole)
D{2} = [1/3 2/3]';   % Exteroceptive states (arousing, non-arousing visual object)
% D{2} = [1 0]';
% Likelihood
%--------------------------------------------------------------------------
A{1}(:,:,1) = [zE*zI         zE*zI         0  ; % Diastole + arousing visual stimulus
               zE*(1-zI)     zE*(1-zI)     1/2; % Systole  + arousing visual stimulus
               zI*(1-zE)     zI*(1-zE)     0  ; % Diastole + non-arousing visual stimulus
               (1-zE)*(1-zI) (1-zE)*(1-zI) 1/2];% Systole  + non-arousing visual stimulus
               
A{1}(:,:,2) = [zI*(1-zE)     zI*(1-zE)     0   ;% Diastole + arousing visual stimulus
               (1-zE)*(1-zI) (1-zE)*(1-zI) 1/2 ;% Systole  + arousing visual stimulus
               zE*zI         zE*zI         0   ;% Diastole + non-arousing visual stimulus
               zE*(1-zI)     zE*(1-zI)     1/2];% Systole  + non-arousing visual stimulus

% Transition probabilities
%--------------------------------------------------------------------------
% Autonomic control states
B{1}(:,:,1) = [(1-wI)/2 (1-wI)/2 wI      ;  % Parasympathetic
               wI       (1-wI)/2 (1-wI)/2;
               (1-wI)/2  wI      (1-wI)/2];
B{1}(:,:,2) = [(1-wI)/2 (1-wI)/2 wI      ;
               (1-wI)/2 (1-wI)/2 (1-wI)/2;
               wI       wI       (1-wI)/2]; % Sympathetic
% Visual object
B{2} = eye(2)*wE+(ones(2)-eye(2))*(1-wE);
           
% Prior preferences
%--------------------------------------------------------------------------
C{1} = [0 xi xi 0]';

% MDP
%--------------------------------------------------------------------------
MDP.A = A;
MDP.B = B;
MDP.C = C;
MDP.D = D;
MDP.T = 40;
% MDP.E = [11/20 9/20]; % prior belief that parasympathetic policy more probable

mdp = spm_MDP_check(MDP);
MDP = spm_MDP_VB_X(mdp);

spm_figure('GetWin','Figure 1');
spm_MDP_VB_trial(MDP)

spm_figure('GetWin','Figure 2');
MDP_ECG_Animation(MDP)



function MDP_ECG_Animation(MDP)

A = Spider_Flower;
x = [];
for j = 1:size(MDP.xn{2},4)
    x(1,end+1:end+size(MDP.xn{2},1)) = MDP.xn{2}(:,1,j,j);
end

ecg   = [0.3*sin(0:0.5:pi) zeros(1,10) -0.5 0:3:5 7:-3:-3 zeros(1,10)];
trace = [];
H = MDP.s(1,:)==3;


for i = 1:length(H)
    if H(i)
        trace(end+1:end+length(ecg)) = ecg;
    elseif i > 1
        if H(i-1)
            trace(end+1:end+length(ecg)) = [0.8*sin(0:0.1:pi) zeros(1,length(ecg)-32)];
        else
            trace(end+1:end+length(ecg)) = zeros(1,length(ecg));
        end
    else
        trace(end+1:end+length(ecg)) = zeros(1,length(ecg));
    end
end

c     = length(x)/length(trace);
for i = 1:length(trace)
    subplot(2,2,1)
    if i < 151
        plot(trace(1:i),'r','LineWidth',3)
    else
        plot(trace(i-150:i),'r','LineWidth',3)
    end
    axis square, axis off
    axis([0 150 -3 8])
    
    subplot(2,2,2)
    ind = max(round(i*c),1);
    I = A.s*x(ind)+A.f*(1-x(ind));
    imagesc(I)
    axis square, axis off
    pause(0.001)
end


function A = Spider_Flower
pth = fileparts(mfilename('fullpath'));
A.s = imread(fullfile(pth,'spider.png'));
A.f = imread(fullfile(pth,'flower.png'));

A.f(end+1,:) = 225;
A.f(:,end+1:end+10,:) = 255;
A.f = [255*ones(size(A.f,1),10,size(A.f,3)) A.f];
