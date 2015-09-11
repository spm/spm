function [] = mci_plot_outputs (M,G)
% Plot outputs
% FORMAT [] = mci_plot_outputs (M,G)
% 
% M     Model
% G     Data
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: mci_plot_outputs.m 6548 2015-09-11 12:39:47Z will $

lw=2;
h=figure;
set(h,'Name','Outputs');
plot(M.t,G,'LineWidth',lw);
for i=1:M.l,
    outstr{i}=sprintf('y(%d)',i);
end
legend(outstr);
xlabel('Time');
grid on