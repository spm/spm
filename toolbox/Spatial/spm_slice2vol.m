function out = spm_slice2vol(job)
% Slice-to-volume alignment job
%__________________________________________________________________________
% Copyright (C) 2021 Wellcome Centre for Human Neuroimaging

% John Ashburner
% $Id: spm_slice2vol.m 8055 2021-02-08 14:28:18Z john $

% Define the output filenames
if job.fwhm>0, prefix = 'sr'; else, prefix = 'r'; end
rfiles = cell(numel(job.images),1);
for i=1:numel(rfiles)
    rfiles{i} = spm_file(job.images{i}, 'prefix',prefix);
end
out.rfiles   = rfiles;
out.rmean{1} = spm_file(job.images{1}, 'prefix','mean', 'number','');
out.rparams  = spm_file(job.images{1}, 'prefix','rp_', 'ext','.mat');

% Estimate the motion
[Q,mu,Mmu,slice_o] = spm_slice2vol_estimate(strvcat(job.images),...
    job.slice_code, job.sd, job.sd*pi/180);

% Save parameters
sd     = job.sd;
images = strvcat(job.images);
save(out.rparams, 'Q','slice_o','sd','images');


plot_parameters(Q,slice_o);


% Write out the mean
dat = file_array(out.rmean{1},size(mu),'int16',0,max(mu(:))/32767,0.0);
nii = nifti;
nii.dat         = dat;
nii.mat         = Mmu;
nii.mat0        = Mmu;
nii.mat_intent  = 'Aligned';
nii.mat0_intent = 'Aligned';
nii.descrip     = 'slice2vol: mean';
create(nii);
nii.dat(:,:,:)  = mu;

% Reslice the images
spm_slice2vol_reslice(strvcat(job.images),Q, job.fwhm);
%==========================================================================

%==========================================================================
function plot_parameters(Q,slice_o)
fg = spm_figure('FindWin','Graphics');
if isempty(fg), return; end

%-Display results: translation and rotation over time series
%--------------------------------------------------------------------------
spm_figure('Clear','Graphics');
ax = axes('Position',[0.1 0.8 0.9 0.1],'Parent',fg,'Visible','off');
set(get(ax,'Title'),'String','Slice to volume alignment',...
    'FontSize',16,'FontWeight','Bold','Visible','on');

% Translations, rotations and time points
Qt = reshape(Q(1:3,slice_o,:),[3 size(Q,2)*size(Q,3)])';
Qr = reshape(Q(4:6,slice_o,:),[3 size(Q,2)*size(Q,3)])'*180/pi;
t  = ((1:size(Qt,1))'-1)/size(Q,2);


ax = axes('Position',[0.1 0.5 0.8 0.3],'Parent',fg,'XGrid','on','YGrid','on',...
    'NextPlot','replacechildren','ColorOrder',[0 0 1;0 0.5 0;1 0 0]);
plot(t,Qt,'-','Parent',ax)
s  = {'x translation','y translation','z translation'};
legend(ax, s, 'Location','Best')
set(get(ax,'Title'),'String','translation','FontSize',16,'FontWeight','Bold');
set(get(ax,'Xlabel'),'String','image');
set(get(ax,'Ylabel'),'String','mm');

ax = axes('Position',[0.1 0.1 0.8 0.3],'Parent',fg,'XGrid','on','YGrid','on',...
    'NextPlot','replacechildren','ColorOrder',[0 0 1;0 0.5 0;1 0 0]);
plot(t,Qr,'-','Parent',ax)
s  = {'pitch','roll','yaw'};
legend(ax, s, 'Location','Best')
set(get(ax,'Title'),'String','rotation','FontSize',16,'FontWeight','Bold');
set(get(ax,'Xlabel'),'String','image');
set(get(ax,'Ylabel'),'String','degrees');

zoom(fg,'on');
drawnow

%-Print realigment parameters
%--------------------------------------------------------------------------
spm_print;



