function spm_fnirs_viewer_sensor(R)
% Display optode/channel positions on the rendered brain surface
% 
% FORMAT spm_fnirs_viewer_sensor(R)
%
% R - structure array containing optode/channel positions 
%    - This structure can be obtained using the SPM-fNIRS toolbox
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny & Sungho Tak
% $Id: spm_fnirs_viewer_sensor.m 6754 2016-03-25 06:44:58Z will $

view = 2; % default view of brain (dorsal view) 

%------------------------------------------------------------------------------
% display channel and optodes positions on the rendered brain 
figure('name', 'Optode and Channel Positions');
brain = R.rend{view}.ren; 
imagesc(brain, [0 2]); 
colormap gray; 
axis image
axis off 

%------------------------------------------------------------------------------
% channels 
% set font and background colors 
bcolor = [1 1 1]; fcolor = [1 0 0]; % ROI

indx = find(sum(R.ch.xy{view}) ~= 0); 
nch = length(indx);

for i = 1:nch 
    r = R.ch.xy{view}(1,indx(i));
    c = R.ch.xy{view}(2,indx(i)); 
    
    ch = R.ch.label(indx(i)); 
    text(c, r, num2str(ch), 'color', fcolor, 'FontWeight', 'bold', 'FontSize', 7, 'HorizontalAlignment', 'center', 'BackgroundColor', bcolor, 'Margin', 0.5); 
end
    
%------------------------------------------------------------------------------
% sources and detectors 
indx = find(sum(R.s.xy{view}) ~= 0); 
ns = length(indx);
for i = 1:ns
    r = R.s.xy{view}(1,indx(i));
    c = R.s.xy{view}(2,indx(i)); 
    text(c, r, 'o', 'color', 'b', 'FontWeight', 'bold', 'FontSize', 10, 'HorizontalAlignment', 'center'); 
end

indx = find(sum(R.d.xy{view}) ~= 0); 
nd = length(indx);
for i = 1:nd
    r = R.d.xy{view}(1,indx(i));
    c = R.d.xy{view}(2,indx(i)); 
    text(c, r, 'x', 'color', 'g', 'FontWeight', 'bold', 'FontSize', 10,  'HorizontalAlignment', 'center');
end
