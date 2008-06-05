function [ZI,f] = spm_eeg_plotScalpData(Z,pos,ChanLabel,in)
%
% function [ZI,f] = spm_eeg_plotScalpData(Z,pos,ChanLabel,in)
%
% This function creates a figure whose purpose is to display an
% interpolation of the sensor data on the scalp (an image)
% IN:
%   - Z: the data matrix at the sensors
%   - pos: the positions of the sensors
%   - ChanLabel: the names of the sensors
%   - in: a structure containing some informations related to the main
%   PRESELECTDATA window. This entry is not necessary.
% OUT:
%   - ZI: an image of interpolated data onto the scalp.
%   - f: the handle of the figure which displays the interpolated data.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jean Daunizeau
% $Id: spm_eeg_plotScalpData.m 1792 2008-06-05 15:43:54Z jean $

if ~exist('in','var') || isempty(in) == 1
    in = [];
    clim = [min(Z(:))-1,max(Z(:))];
    figName = ['Image Scalp data'];
    deleteFcn = [];
else
    clim = [in.min,in.max];
    if ~isfield(in,'trN')
        figName = ['Image Scalp data: ',in.type,' sensors'];
    else
        figName = ['Image Scalp data: ',in.type,' sensors, trial #',num2str(in.trN),'.'];
    end
    deleteFcn = ['try;delete(get(gcf,''userdata''));end'];
end

if ~isequal(size(pos,2),length(ChanLabel))
    pos = pos';
end
nD = size(pos,1);
if nD ~= 2
    % get 2D positions from 3D positions
   xyz = pos;
   [pos] = get2Dfrom3D(xyz);
   pos = pos';
end

% exclude channels ?
goodChannels = find(~isnan(pos(1,:)));
pos = pos(:,goodChannels);
Z = Z(goodChannels,:);
ChanLabel = ChanLabel(goodChannels);

xmin = min(pos(1,:));
xmax = max(pos(1,:));
dx = (xmax-xmin)./100;
ymin = min(pos(2,:));
ymax = max(pos(2,:));
dy = (ymax-ymin)./100;
x = xmin:dx:xmax;
y = ymin:dy:ymax;
[XI,YI] = meshgrid(x,y);
ZI = griddata(pos(1,:)',pos(2,:)',double(Z'),XI,YI);

f=figure;
set(f,'name',figName,'deleteFcn',deleteFcn);
hi = imagesc(flipud(ZI));
hold on
[C,hc] = contour(flipud(ZI));
set(hc,'linecolor',0.5.*ones(3,1))
hold off
caxis(clim);
col = colormap;
col(1,:) = .8*ones(3,1);
colormap(col)
colorbar
axis off
axis equal
axis tight

fpos = pos;
fpos(1,:) = fpos(1,:) - xmin;
fpos(2,:) = fpos(2,:) - ymin;
fpos(1,:) = fpos(1,:)./(dx);
fpos(2,:) = fpos(2,:)./(dy);
fpos(2,:) = 100-fpos(2,:);  % for display purposes (flipud imagesc)

figure(f);
hold on;
hp = plot(fpos(1,:),fpos(2,:),'ko');
ht = text(fpos(1,:),fpos(2,:),ChanLabel);
set(ht,'visible','off')
axis image

d.interp.XI = XI;
d.interp.YI = YI;
d.interp.pos = pos;
d.f = f;
d.pos = fpos;
d.goodChannels = goodChannels;
d.ChanLabel = ChanLabel;
d.hp = hp;
d.ht = ht;
d.hi = hi;
d.in = in;

d.hsp = uicontrol('style','pushbutton','callback',{@dosp},...
    'BusyAction','cancel',...
    'Interruptible','off',...
    'position',[10    120    80    20],...
    'string','channel pos');
d.hsn = uicontrol('style','pushbutton','callback',{@dosn},...
    'BusyAction','cancel',...
    'Interruptible','off',...
    'position',[10    80    80    20],...
    'string','channel names');
if ~isempty(in)
    ud = get(in.handles.hfig,'userdata');
    nT = ud.Nsamples;
    d.hti = uicontrol('style','text',...
        'string',[num2str(in.gridTime(in.x)),' (',...
        in.unit,')'],...
        'position',[10    10    120    20]);
    d.hts = uicontrol('style','slider',...
        'Position',[130 10 30 20],'min',1,'max',nT,...
        'value',in.x,'sliderstep',[1./(nT-1) 1./(nT-1)],...
        'callback',{@doChangeTime},...
        'BusyAction','cancel',...
        'Interruptible','off');
    set(d.hti,'userdata',d);
    set(d.hts,'userdata',d);
end
set(d.hsp,'userdata',d);
set(d.hsn,'userdata',d);




function dosp(btn,evd)
dbstop if error
d = get(btn,'userdata');
v = get(d.hp,'visible');
switch v
    case 'on'
        set(d.hp,'visible','off');
    case 'off'
        set(d.hp,'visible','on');
end


function dosn(btn,evd)
dbstop if error
d = get(btn,'userdata');
v = get(d.ht(1),'visible');
switch v
    case 'on'
        set(d.ht,'visible','off');
    case 'off'
        set(d.ht,'visible','on');
end


function doChangeTime(btn,evd)
dbstop if error
d = get(btn,'userdata');
v = get(btn,'value');
% get data
D = get(d.in.handles.hfig,'userdata');
if ~isfield(d.in,'trN')
    trN = 1;
else
    trN = d.in.trN;
end
Z = D.data.y(d.in.ind,v,trN);
Z = Z(d.goodChannels);
clear ud;
% interpolate data
ZI = griddata(d.interp.pos(1,:),d.interp.pos(2,:),double(Z),d.interp.XI,d.interp.YI);
% update data display
set(d.hi,'Cdata',flipud(ZI));
% update time index display
set(d.hti,'string',[num2str(d.in.gridTime(v)),' (',...
        d.in.unit,')']);
% update display marker position
try;set(d.in.hl,'xdata',[v;v]);end
hf=findobj(gca,'type','hggroup');
delete(hf)
hold on
[C,hc] = contour(flipud(ZI));
set(hc,'linecolor',[0.5.*ones(3,1)])
hold off
drawnow
axis image


function [xy] = get2Dfrom3D(xyz)
% function [xy] = get2Dfrom3D(xyz)
% This function is used to flatten 3D sensor positions onto the 2D plane
% using a modified spherical projection operation.
% It is used to visualize channel data.
% IN:
%   - xyz: the carthesian sensor position in 3D space
% OUT:
%   - xy: the (x,y) carthesian coordinates of the sensors after projection
%   onto the best-fitting sphere
dbstop if error
if size(xyz,2) ~= 3
    xyz = xyz';
end
% exclude channels ?
badChannels = find(isnan(xyz(:,1)));
goodChannels = find(isnan(xyz(:,1))~=1);
xyz = xyz(goodChannels,:);
% Fit sphere to 3d sensors and center frame
[C,R,out] = fitSphere(xyz(:,1),xyz(:,2),xyz(:,3));
xyz = xyz - repmat(C,size(xyz,1),1);
% apply transformation using spherical coordinates
[TH,PHI,RAD] = cart2sph(xyz(:,1),xyz(:,2),xyz(:,3));
TH = TH - mean(TH);
[X,Y,Z] = sph2cart(TH,zeros(size(TH)),RAD.*(cos(PHI+pi./2)+1));
xy = [X(:),Y(:)];
