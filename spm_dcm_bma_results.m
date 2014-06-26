function spm_dcm_bma_results(BMS,method)
% Plot histograms from BMA for selected modulatory and driving input
% FORMAT spm_dcm_bma_results(BMS,method)
% BMS        - BMS.mat file
% method     - inference method (FFX or RFX)
%__________________________________________________________________________
% Copyright (C) 2009-2014 Wellcome Trust Centre for Neuroimaging

% Maria Joao
% $Id: spm_dcm_bma_results.m 6067 2014-06-26 15:33:30Z guillaume $


%-Load BMS file
%--------------------------------------------------------------------------
if ~nargin
    [fname, sts] = spm_select(1,'^BMS.mat$','select BMS.mat file');
    if ~sts, return; end
else
    fname = BMS;
end

load(fname)

%-Check BMS/BMA method used
%--------------------------------------------------------------------------
if nargin < 2
    ff = fieldnames(BMS.DCM);
    Nff = numel(ff);
    if Nff==2
        method = spm_input('Inference method','+1','b','FFX|RFX',['ffx';'rfx']);
    else % pick the one available if only one method
        method = char(ff);
    end
end

%-Select method
%--------------------------------------------------------------------------
if isfield(BMS.DCM,method)
    switch method
        case 'ffx'
            if isempty(BMS.DCM.ffx.bma)
                error('No BMA analysis for FFX in BMS file.');
            else

                Nsamp = BMS.DCM.ffx.bma.nsamp;
                amat  = BMS.DCM.ffx.bma.a;
                bmat  = BMS.DCM.ffx.bma.b;
                cmat  = BMS.DCM.ffx.bma.c;
                dmat  = BMS.DCM.ffx.bma.d;
            end
            disp('Loading model space...')
            load(BMS.DCM.ffx.data)
            load(subj(1).sess(1).model(1).fname)

        case 'rfx'
            if isempty(BMS.DCM.rfx.bma)
                error('No BMA analysis for RFX in BMS file.');
            else
                Nsamp = BMS.DCM.rfx.bma.nsamp;
                amat  = BMS.DCM.rfx.bma.a;
                bmat  = BMS.DCM.rfx.bma.b;
                cmat  = BMS.DCM.rfx.bma.c;
                dmat  = BMS.DCM.rfx.bma.d;
            end
            disp('Loading model space...')
            load(BMS.DCM.rfx.data)
            load(subj(1).sess(1).model(1).fname)
    end
else
    error('No %s analysis in current BMS.mat file.',method);
end

%-Number of regions, mod. inputs and names
%--------------------------------------------------------------------------
n  = size(amat,2); % #region
m  = size(bmat,3); % #drv/mod inputs
mi = size(cmat,2);

% Look for modulatory inputs
mod_input = [];
for ii=1:m
    % look for bits of B not full of zeros
    tmp = squeeze(bmat(:,:,ii,:));
    if any(tmp(:))
        mod_input = [mod_input ii];
    end
end
% Look for effective driving inputs
drive_input = [];
for ii=1:m
    % look for bits of not full of zeros
    tmp = any(cmat(:,ii,:));
    if sum(tmp)
        drive_input = [drive_input ii];
    end
end

% Non linear model ? If so find the driving regions
if ~isempty(dmat)
    nonLin = 1;
    mod_reg = [];
    for ii=1:n
        % look for bits of D not full of zeros
        tmp = squeeze(dmat(:,:,ii,:));
        if any(tmp(:))
            mod_reg = [mod_reg ii];
        end
    end
else
    nonLin = 0;
    mod_reg = [];
end

if isfield(DCM.Y,'name')
    for i=1:n
        region(i).name = DCM.Y.name{i};
    end
else
    for i=1:n
        str            = sprintf('Region %d',i);
        region(i).name = spm_input(['Name for ',str],'+1','s');
    end
end

bins   = Nsamp/100;

%-Intrinsic connection density
%--------------------------------------------------------------------------
F  = spm_figure('GetWin','Graphics');
set(F,'name',sprintf('%s: %s',spm('Version'),'BMA Results'));
FS = spm('FontSizes');

usd.amat        = amat;
usd.bmat        = bmat;
usd.cmat        = cmat;
usd.dmat        = dmat;

usd.region      = region;
usd.n           = n;
usd.m           = m;
usd.ni          = mi;
usd.FS          = FS;
usd.drive_input = drive_input;
usd.mod_input   = mod_input;
if nonLin
    usd.mod_reg = mod_reg;
end
usd.bins        = bins;
usd.Nsamp       = Nsamp;

set(F,'userdata',usd);
clf(F);

labels = {'a: int.'};
callbacks = {@plot_a};
for ii = mod_input
    labels{end+1} = ['b: mod. i#',num2str(ii)];
    callbacks{end+1} = @plot_b;
end
for ii = drive_input
    labels{end+1} = ['c: drv. i#',num2str(ii)];
    callbacks{end+1} = @plot_c;
end

if nonLin
    for ii = mod_reg
        labels{end+1} = ['d: mod. r#',num2str(ii)];
        callbacks{end+1} = @plot_d;
    end
end   

[handles] = spm_uitab(F,labels,callbacks,'BMA_parameters',1);

set(handles.htab,'backgroundcolor',[1 1 1])
set(handles.hh,'backgroundcolor',[1 1 1])
set(handles.hp,'HighlightColor',0.8*[1 1 1])
set(handles.hp,'backgroundcolor',[1 1 1])

plot_a(F);


%==========================================================================
function plot_a(F)

if ~nargin
    F = get(gco,'parent');
end

H = findobj(F,'tag','BMA_parameters','type','uipanel');

hc = intersect(findobj('tag','bma_results'),get(H,'children'));
if ~isempty(hc)
    delete(hc)
end

ud = get(F,'userdata');

titlewin = 'BMA: intrinsic connections (a)';
hTitAx = axes('Parent',H,'Position',[0.2,0.04,0.6,0.02],...
    'Visible','off','tag','bma_results');
text(0.55,0,titlewin,'Parent',hTitAx,'HorizontalAlignment','center',...
    'VerticalAlignment','baseline','FontWeight','Bold','FontSize',ud.FS(12))

for i=1:ud.n
    for j=1:ud.n
        k=(i-1)*ud.n+j;
        subplot(ud.n,ud.n,k);
        if (i==j)
            axis off
        else
            hist(squeeze(ud.amat(i,j,:)),ud.bins,'r');
            amax = max(abs(ud.amat(i,j,:)))*1.05; % enlarge scale by 5%
            if amax > 0
                xlim([-amax amax])
            else % case where parameter is constrained to be 0.
                xlim([-10 10])
            end
            set(gca,'YTickLabel',[]);
            set(gca,'FontSize',12);
            title(sprintf('%s to %s',ud.region(j).name,ud.region(i).name));
        end
    end
end

%==========================================================================
function plot_b

hf = get(gco,'parent');
ud = get(hf,'userdata');
H  = findobj(hf,'tag','BMA_parameters','type','uipanel');

hc = intersect(findobj('tag','bma_results'),get(H,'children'));
if ~isempty(hc)
    delete(hc)
end

% spot the bmod input index from the fig name
ht = intersect(findobj('style','pushbutton'),get(hf,'children'));
ht = findobj(ht,'flat','Fontweight','bold');
t_str = get(ht,'string');
b_ind = str2num(t_str(strfind(t_str,'#')+1:end));
i_mod = find(ud.mod_input==b_ind);

titlewin = ['BMA: modulatory connections (b',num2str(b_ind),')'];
hTitAx = axes('Parent',H,'Position',[0.2,0.04,0.6,0.02],...
    'Visible','off','tag','bma_results');
text(0.55,0,titlewin,'Parent',hTitAx,'HorizontalAlignment','center',...
    'VerticalAlignment','baseline','FontWeight','Bold','FontSize',ud.FS(12))

for i=1:ud.n
    for j=1:ud.n
        k=(i-1)*ud.n+j;
        subplot(ud.n,ud.n,k);
        if (i==j)
            axis off
        else
            hist(squeeze(ud.bmat(i,j,ud.mod_input(i_mod),:)),ud.bins,'r');
            bmax = max(abs(ud.bmat(i,j,ud.mod_input(i_mod),:)))*1.05; % enlarge scale by 5%
            if bmax > 0
                xlim([-bmax bmax])
            else % case where parameter is constrained to be 0.
                xlim([-10 10])
            end
            set(gca,'YTickLabel',[]);
            set(gca,'FontSize',12);
            title(sprintf('%s to %s',ud.region(j).name,ud.region(i).name));
        end
    end
end

%==========================================================================
function plot_c

hf = get(gco,'parent');
ud = get(hf,'userdata');
H  = findobj(hf,'tag','BMA_parameters','type','uipanel');

hc = intersect(findobj('tag','bma_results'),get(H,'children'));
if ~isempty(hc)
    delete(hc)
end

% spot the c_drv input index from the fig name
ht = intersect(findobj('style','pushbutton'),get(hf,'children'));
ht = findobj(ht,'flat','Fontweight','bold');
t_str = get(ht,'string');
c_ind = str2num(t_str(strfind(t_str,'#')+1:end));
i_drv = find(ud.drive_input==c_ind);

titlewin = ['BMA: input connections (c',num2str(c_ind),')'];
hTitAx = axes('Parent',H,'Position',[0.2,0.04,0.6,0.02],...
    'Visible','off','tag','bma_results');
text(0.55,0,titlewin,'Parent',hTitAx,'HorizontalAlignment','center',...
    'VerticalAlignment','baseline','FontWeight','Bold','FontSize',ud.FS(12))

for j=1:ud.n
    subplot(1,ud.n,j);
    if length(find(ud.cmat(j,ud.drive_input(i_drv),:)==0))==ud.Nsamp
        plot([0 0],[0 1],'k');
    else
        hist(squeeze(ud.cmat(j,ud.drive_input(i_drv),:)),ud.bins,'r');
        cmax = max(abs(ud.cmat(j,ud.drive_input(i_drv),:)))*1.05; % enlarge scale by 5%
        if cmax > 0
            xlim([-cmax cmax])
        else % case where parameter is constrained to be 0.
            xlim([-10 10])
        end
    end
    set(gca,'YTickLabel',[]);
    set(gca,'FontSize',12);
    title(sprintf('%s ',ud.region(j).name));
end

%==========================================================================
function plot_d

hf = get(gco,'parent');
ud = get(hf,'userdata');
H  = findobj(hf,'tag','BMA_parameters','type','uipanel');

hc = intersect(findobj('tag','bma_results'),get(H,'children'));
if ~isempty(hc)
    delete(hc)
end

% spot the d_reg input index from the fig name
ht = intersect(findobj('style','pushbutton'),get(hf,'children'));
ht = findobj(ht,'flat','Fontweight','bold');
t_str = get(ht,'string');
d_ind = str2num(t_str(strfind(t_str,'#')+1:end));
i_mreg = find(ud.mod_reg==d_ind);

titlewin = ['BMA: non-linear connections (d',num2str(d_ind),')'];
hTitAx = axes('Parent',H,'Position',[0.2,0.04,0.6,0.02],...
    'Visible','off','tag','bma_results');
text(0.55,0,titlewin,'Parent',hTitAx,'HorizontalAlignment','center',...
    'VerticalAlignment','baseline','FontWeight','Bold','FontSize',ud.FS(12))

for i=1:ud.n
    for j=1:ud.n
        k=(i-1)*ud.n+j;
        subplot(ud.n,ud.n,k);
        if (i==j)
            axis off
        else
            hist(squeeze(ud.dmat(i,j,ud.mod_reg(i_mreg),:)),ud.bins,'r');
            dmax = max(abs(ud.dmat(i,j,ud.mod_reg(i_mreg),:)))*1.05; % enlarge scale by 5%
            if dmax > 0
                xlim([-dmax dmax])
            else % case where parameter is constrained to be 0.
                xlim([-10 10])
            end
            set(gca,'YTickLabel',[]);
            set(gca,'FontSize',12);
            title(sprintf('%s to %s',ud.region(j).name,ud.region(i).name));
        end
    end
end
