function spm_dcm_bma_results (BMS,mod_in,drive_in,method)
% Plot histograms from BMA for selected modulatory and driving input
% FORMAT spm_dcm_bma_results (BMS,mod_in,drive_in,method)
%
% Input:
% BMS        - BMS.mat file 
% mod_in     - modulatory input
% drive_in   - driving input
% method     - inference method (FFX or RFX)
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Maria Joao
% $Id: spm_dcm_bma_results.m 3632 2009-12-11 09:58:31Z maria $

if nargin < 6
    % function called without parameters (e.g. via GUI)
    %----------------------------------------------------------------------
    Finter = spm_figure('GetWin','Interactive');
    spm_clf(Finter);
    set(Finter,'name','Dynamic Causal Modeling');
    
    header      = get(Finter,'Name');
    WS          = spm('WinScale');
    fname       = spm_select([1 1],'^BMS.mat$','select BMS.mat file');
    
    mod_input   = spm_input('Select modulatory input ? ',1,'r',[],1);
    drive_input = spm_input('Select driving input ? ','+1','r',[],1);
    method      = spm_input('Inference method','+1','b','FFX|RFX',['ffx';'rfx']);    
else
    % use function arguments
    %----------------------------------------------------------------------
    mod_input   = mod_in;
    drive_input = drive_in;
    fname       = BMS;   
end

% load BMS file
%--------------------------------------------------------------------------
load(fname)

nonLin = 0;
% select method
%--------------------------------------------------------------------------
if isfield(BMS.DCM,method)
    switch method
        case 'ffx'
                if isempty(BMS.DCM.ffx.bma)
                    error('No BMA analysis for FFX in BMS file!');
                else
                    Nsamp = BMS.DCM.ffx.bma.nsamp;
                    amat  = BMS.DCM.ffx.bma.a;
                    bmat  = BMS.DCM.ffx.bma.b;
                    cmat  = BMS.DCM.ffx.bma.c;
                    if isfield(BMS.DCM.ffx.bma,'d')
                        dmat = BMS.DCM.ffx.bma.d;
                        mod_reg = spm_input('Select modulating region ? ','+1','r',[],1);
                        nonLin = 1;
                    end
                    
                end
                disp('Loading model space...')
                load(BMS.DCM.ffx.data)
                load(subj(1).sess(1).model(1).fname)
            
        case 'rfx'
                if isempty(BMS.DCM.rfx.bma)
                    error('No BMA analysis for RFX in BMS file!');
                else
                    Nsamp = BMS.DCM.rfx.bma.nsamp;
                    amat = BMS.DCM.rfx.bma.a;
                    bmat = BMS.DCM.rfx.bma.b;
                    cmat = BMS.DCM.rfx.bma.c;
                    if isfield(BMS.DCM.rfx.bma,'d')
                       dmat = BMS.DCM.rfx.bma.d;
                       mod_reg   = spm_input('Select modulating region ? ','+1','r',[],1);
                       nonLin = 1;
                    end
                    
                end
                disp('Loading model space...')
                load(BMS.DCM.rfx.data)
                load(subj(1).sess(1).model(1).fname)
    end
else
    msgbox(sprintf('No %s analysis in current BMS.mat file!',method))
    return
end

% number of regions, mod. inputs and names
%--------------------------------------------------------------------------
n  = size(amat,2);
m  = size(bmat,3);
mi = size(cmat,2);

% check if input is correct
if mod_input > m || drive_input > mi
    error('Incorrect choice for driving or modulatory input!');
end

if isfield(DCM.Y,'name')
    for i=1:n,
        region(i).name = DCM.Y.name{i};
    end
else
    for i=1:n,
        str            = sprintf('Region %d',i);
        region(i).name = spm_input(['Name for ',str],'+1','s');
    end
end

bins = Nsamp/100;

% intrinsic connection density
%--------------------------------------------------------------------------
F  = spm_figure('GetWin','Graphics','BMA: results');
FS = spm('FontSizes');

usd.amat        = amat;
usd.bmat        = bmat;
usd.cmat        = cmat;
if nonLin
    usd.dmat    = dmat;
    usd.mod_reg = drive_input;
end
usd.region      = region;
usd.n           = n;
usd.m           = m;
usd.ni          = mi;
usd.FS          = FS;
usd.drive_input = drive_input;
usd.mod_input   = mod_input;
usd.bins        = bins;
usd.Nsamp       = Nsamp;

set(F,'userdata',usd);
clf(F);

if nonLin
    labels = {'a: int.','b: mod.','c: inp.', 'd: non.'};
    callbacks = {@plot_a,@plot_b,@plot_c,@plot_d};
else
    labels = {'a: int.','b: mod.','c: inp.'};
    callbacks = {@plot_a,@plot_b,@plot_c};
end

[handles] = spm_uitab(F,labels,callbacks,'BMA_parameters',1);

set(handles.htab,'backgroundcolor',[1 1 1])
set(handles.hh,'backgroundcolor',[1 1 1])
set(handles.hp,'HighlightColor',0.8*[1 1 1])
set(handles.hp,'backgroundcolor',[1 1 1])

feval(@plot_a,F)

end

function plot_a(F)

try
    F;
catch
    F = get(gco,'parent');
end

hc = intersect(findobj('tag','bma_results'),get(F,'children'));
if ~isempty(hc)
    delete(hc)
end

ud = get(F,'userdata');

titlewin = 'BMA: intrinsic connections (a)';
hTitAx = axes('Parent',F,'Position',[0.2,0.04,0.6,0.02],...
    'Visible','off','tag','bma_results');
text(0.55,0,titlewin,'Parent',hTitAx,'HorizontalAlignment','center',...
    'VerticalAlignment','baseline','FontWeight','Bold','FontSize',ud.FS(12))

for i=1:ud.n,
    for j=1:ud.n,
        k=(i-1)*ud.n+j;
        subplot(ud.n,ud.n,k);
        if (i==j)
            axis off
        else
            hist(ud.amat(i,j,:),ud.bins,'r');
            amax = max(ud.amat(i,j,:));
            if amax > 0
                xlim([-amax amax])
            else
                if amax < 0
                   amin = min(ud.amat(i,j,:));
                   xlim([amin amax])
                else
                    xlim([-10 10])
                end
            end
            set(gca,'YTickLabel',[]);
            set(gca,'FontSize',12);
            title(sprintf('%s to %s',ud.region(j).name,ud.region(i).name));
        end
    end
end

end

function plot_b()

hf = get(gco,'parent');
ud = get(hf,'userdata');

hc = intersect(findobj('tag','bma_results'),get(hf,'children'));
if ~isempty(hc)
    delete(hc)
end

titlewin = 'BMA: modulatory connections (b)';
hTitAx = axes('Parent',hf,'Position',[0.2,0.04,0.6,0.02],...
    'Visible','off','tag','bma_results');
text(0.55,0,titlewin,'Parent',hTitAx,'HorizontalAlignment','center',...
    'VerticalAlignment','baseline','FontWeight','Bold','FontSize',ud.FS(12))

for i=1:ud.n,
    for j=1:ud.n,
        k=(i-1)*ud.n+j;
        subplot(ud.n,ud.n,k);
        if (i==j)
            axis off
        else
            hist(ud.bmat(i,j,ud.mod_input,:),ud.bins,'r');
            bmax = max(ud.bmat(i,j,ud.mod_input,:));
            if bmax > 0
                xlim([-bmax bmax])
            else
                if bmax < 0
                   bmin = min(ud.bmat(i,j,ud.mod_input,:));
                   xlim([bmin bmax])
                else
                    xlim([-10 10])
                end
            end
            set(gca,'YTickLabel',[]);
            set(gca,'FontSize',12);
            title(sprintf('%s to %s',ud.region(j).name,ud.region(i).name));
        end
    end
end


end

function plot_c()

hf = get(gco,'parent');
ud = get(hf,'userdata');

hc = intersect(findobj('tag','bma_results'),get(hf,'children'));
if ~isempty(hc)
    delete(hc)
end

titlewin = 'BMA: input connections (c)';
hTitAx = axes('Parent',hf,'Position',[0.2,0.04,0.6,0.02],...
    'Visible','off','tag','bma_results');
text(0.55,0,titlewin,'Parent',hTitAx,'HorizontalAlignment','center',...
    'VerticalAlignment','baseline','FontWeight','Bold','FontSize',ud.FS(12))

for j=1:ud.n,
    subplot(1,ud.n,j);
    if length(find(ud.cmat(j,ud.drive_input,:)==0))==ud.Nsamp
        plot([0 0],[0 1],'k');
    else
        hist(ud.cmat(j,ud.drive_input,:),ud.bins,'r');
        cmax = max(ud.cmat(j,ud.drive_input,:));
        if cmax > 0
           xlim([-cmax cmax])
        else
            if cmax < 0
                cmin = min(ud.cmat(j,ud.drive_input,:));
                xlim([cmin cmax])
            else
                xlim([-10 10])
            end
        end
    end
    set(gca,'YTickLabel',[]);
    set(gca,'FontSize',12);
    title(sprintf('%s ',ud.region(j).name));
end

end

function plot_d()

hf = get(gco,'parent');
ud = get(hf,'userdata');

hc = intersect(findobj('tag','bma_results'),get(hf,'children'));
if ~isempty(hc)
    delete(hc)
end

titlewin = 'BMA: non-linear connections (d)';
hTitAx = axes('Parent',hf,'Position',[0.2,0.04,0.6,0.02],...
    'Visible','off','tag','bma_results');
text(0.55,0,titlewin,'Parent',hTitAx,'HorizontalAlignment','center',...
    'VerticalAlignment','baseline','FontWeight','Bold','FontSize',ud.FS(12))

for i=1:ud.n,
    for j=1:ud.n,
        k=(i-1)*ud.n+j;
        subplot(ud.n,ud.n,k);
        if (i==j)
            axis off
        else
            hist(ud.bmat(i,j,ud.mod_input,:),ud.bins,'r');
            bmax = max(ud.dmat(i,j,ud.mod_input,:));
            if bmax > 0
                xlim([-bmax bmax])
            else
                if bmax < 0
                   bmin = min(ud.dmat(i,j,ud.mod_reg,:));
                   xlim([bmin bmax])
                else
                    xlim([-10 10])
                end
            end
            set(gca,'YTickLabel',[]);
            set(gca,'FontSize',12);
            title(sprintf('%s to %s',ud.region(j).name,ud.region(i).name));
        end
    end
end

end