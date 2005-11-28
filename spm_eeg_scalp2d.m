function [x,y,z] = spm_eeg_scalp2d(D, d)
% interpolation of data d on scalp. Electrode positions are taken from
% channel template file ctf.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% James Kilner & Stefan Kiebel
% $Id: spm_eeg_scalp2d.m 317 2005-11-28 18:31:24Z stefan $

load(fullfile(spm('dir'), 'EEGtemplates', D.channels.ctf));

Cpos = Cpos(:, D.channels.order(D.channels.eeg));

x = min(Cpos(1,:)):0.005:max(Cpos(1,:));
y = min(Cpos(2,:)):0.005:max(Cpos(2,:));

[x1,y1] = meshgrid(x,y);
xp = Cpos(1,:)';
yp = Cpos(2,:)';

z = griddata(xp, yp, d, x1, y1);

figure
surface(x,y,z);
shading('interp')
hold on
plot3(xp, yp, d, 'k.');
