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
% $Id: spm_dcm_bma_results.m 3569 2009-11-13 15:51:07Z guillaume $

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
                    
                end
                load(BMS.DCM.ffx.data(1).sess(1).model(1).fname)
            
        case 'rfx'
                if isempty(BMS.DCM.rfx.bma)
                    error('No BMA analysis for RFX in BMS file!');
                else
                    Nsamp = BMS.DCM.rfx.bma.nsamp;
                    amat = BMS.DCM.rfx.bma.a;
                    bmat = BMS.DCM.rfx.bma.b;
                    cmat = BMS.DCM.rfx.bma.c;
                    
                end
                load(BMS.DCM.rfx.data(1).sess(1).model(1).fname)
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

figure(F);

titlewin = 'BMA: intrinsic connections (a)';
hTitAx = axes('Parent',F,'Position',[0.02 0.97 0.86 0.02],...
        'Visible','off');
text(0.55,0,titlewin,'Parent',hTitAx,'HorizontalAlignment','center',...
        'VerticalAlignment','baseline','FontWeight','Bold','FontSize',FS(12))

for i=1:n,
    for j=1:n,
        k=(i-1)*n+j;
        subplot(n,n,k);
        if (i==j)
            axis off
        else
            hist(amat(i,j,:),bins,'r');
            xlim([-max(amat(i,j,:)) max(amat(i,j,:))])
            set(gca,'YTickLabel',[]);
            set(gca,'FontSize',12);
            title(sprintf('%s to %s',region(j).name,region(i).name));
        end
    end
end

% modulatory connections density
%--------------------------------------------------------------------------
F = spm_figure('CreateWin','Graphics','BMA: results');

figure(F);

titlewin = 'BMA: modulatory connections (b)';
hTitAx = axes('Parent',F,'Position',[0.02 0.97 0.86 0.02],...
        'Visible','off');
text(0.55,0,titlewin,'Parent',hTitAx,'HorizontalAlignment','center',...
        'VerticalAlignment','baseline','FontWeight','Bold','FontSize',FS(12))

for i=1:n,
    for j=1:n,
        k=(i-1)*n+j;
        subplot(n,n,k);
        if (i==j)
            axis off
        else
            hist(bmat(i,j,mod_input,:),bins,'r');
            xlim([-max(bmat(i,j,mod_input,:)) max(bmat(i,j,mod_input,:))])
            set(gca,'YTickLabel',[]);
            set(gca,'FontSize',12);
            title(sprintf('%s to %s',region(j).name,region(i).name));
        end
    end
end

% input connection density
%--------------------------------------------------------------------------
F = spm_figure('CreateWin','Graphics','BMA: results');

titlewin = 'BMA: input connections (c)';
hTitAx = axes('Parent',F,'Position',[0.02 0.97 0.86 0.02],...
        'Visible','off');
text(0.55,0,titlewin,'Parent',hTitAx,'HorizontalAlignment','center',...
        'VerticalAlignment','baseline','FontWeight','Bold','FontSize',FS(12))

for i=1:n,
    subplot(1,n,i);
    if length(find(cmat(i,drive_input,:)==0))==Nsamp
        plot([0 0],[0 1],'k');
    else
        hist(cmat(i,drive_input,:),bins,'r');
        xlim([-max(cmat(i,drive_input,:)) max(cmat(i,drive_input,:))])
    end
    set(gca,'YTickLabel',[]);
    set(gca,'FontSize',12);
    title(sprintf('%s ',region(i).name));
end
