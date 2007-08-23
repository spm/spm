function spm_DFP_plot(QU,pU)
% plots particles for spm_DFP
% FORMAT spm_DFP_plot(QU,Nt)
% FORMAT spm_DFP_plot(QU,pU)
%--------------------------------------------------------------------------
% QU{t}(p).x{d}  - ensemble of hidden states
% QU{t}(p).v{d}  - ensemble of causal states
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id$

% defaults for plotting
%--------------------------------------------------------------------------
clf
try
    pX = pU.x{1};
    pV = pU.v{2};
    Nt = length(pX);
catch
    try Nt = pU; catch, Nt = length(QU); end
end

% time-series specification
%--------------------------------------------------------------------------
nx    = size(QU{1}(1).x{1},1);
nv    = size(QU{1}(1).v{1},1);
np    = length(QU{1});
nt    = length(QU);

% unpack states
%--------------------------------------------------------------------------
for t = 1:nt
    for j = 1:np
        for i = 1:nx
            X{i}(t,j) = QU{t}(j).x{1}(i);
        end
        for i = 1:nv
            V{i}(t,j) = QU{t}(j).v{1}(i);
        end
    end
end

% causes
%--------------------------------------------------------------------------
if nt < 2, return, end

subplot(2,1,1)
for i = 1:nv
    plot(1:nt,V{i},':','Color',[1 1 1]/(2 + i - 1))
    hold on
    plot(1:nt,mean(V{i},2),'--b','LineWidth',2)
    hold on
end
try
    hold on
    plot([1:nt] - 1,pV,'r')
end
hold off
title(sprintf('causal states - %i',nv));
xlabel('time {bins}')
ylabel('states (a.u.)')
grid on
axis square
set(gca,'XLim',[1 Nt])


% hidden states
%--------------------------------------------------------------------------
if ~nx, drawnow, return, end

subplot(2,1,2)
for i = 1:nx
    plot(1:nt,X{i},':','Color',[1 1 1]/(2 + i - 1))
    hold on
    plot(1:nt,mean(X{i},2),'--b','LineWidth',2)
    hold on
end
try
    hold on
    plot([1:nt] - 1,pX,'r')
end
hold off
title(sprintf('hidden states - %i',nx));
xlabel('time {bins}')
ylabel('states (a.u.)')
grid on
axis square
set(gca,'XLim',[1 Nt])
drawnow