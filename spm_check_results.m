function spm_check_results(SPMs,xSPM)
% Display several MIPs in the same figure
% FORMAT spm_check_results(SPMs,xSPM)
% SPMs    - char or cell array of paths to SPM.mat[s]
% xSPM    - structure containing thresholding details, see spm_getSPM.m
%
% Beware: syntax and features of this function are likely to change.
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_check_results.m 4660 2012-02-20 17:39:29Z guillaume $


cwd = pwd;

%-Get input parameter SPMs
%--------------------------------------------------------------------------
if ~nargin || isempty(SPMs)
    [SPMs, sts] = spm_select(Inf,'^SPM\.mat$','Select SPM.mat[s]');
    if ~sts, return; end
end
SPMs = cellstr(SPMs);

%-Get input parameter xSPM
%--------------------------------------------------------------------------
xSPM.swd = spm_file(SPMs{1},'fpath');
try, [xSPM.thresDesc, xSPM.u] = convert_desc(xSPM.thresDesc); end
[tmp, xSPM] = spm_getSPM(xSPM);
if ~isfield(xSPM,'units'), xSPM.units = {'mm' 'mm' 'mm'}; end
[xSPM.thresDesc, xSPM.u] = convert_desc(xSPM.thresDesc);

%-
%--------------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');
spm_figure('Clear','Graphics');

mn = numel(SPMs);
n  = round(mn^0.4);
m  = ceil(mn/n);
w  = 1/n;
h  = 1/m;
ds = (w+h)*0.02;

for ij=1:numel(SPMs)
    i = 1-h*(floor((ij-1)/n)+1);
    j = w*rem(ij-1,n);
    
    xSPM.swd = spm_file(SPMs{ij},'fpath');
    [tmp, myxSPM] = spm_getSPM(xSPM);
    
    hMIPax(ij) = axes('Parent',Fgraph,'Position',[j+ds/2 i+ds/2 w-ds h-ds],'Visible','off');
    hMIPax(ij) = spm_mip_ui(myxSPM.Z,myxSPM.XYZmm,myxSPM.M,myxSPM.DIM,hMIPax(ij),myxSPM.units);
    axis(hMIPax(ij),'image');
    set(findobj(hMIPax(ij),'type','image'),'UIContextMenu',[]);
    hReg = get(hMIPax(ij),'UserData');
    set(hReg.hXr,'visible','off');
    set(hReg.hMIPxyz,'visible','off');
end
linkaxes(hMIPax);

cd(cwd);

%==========================================================================
% function [str, u] = convert_desc(str)
%==========================================================================
function [str, u] = convert_desc(str)
td = regexp(str,'p\D?(?<u>[\.\d]+) \((?<thresDesc>\S+)\)','names');
if isempty(td)
    td = regexp(str,'\w=(?<u>[\.\d]+)','names');
    td.thresDesc = 'none';
end
if strcmp(td.thresDesc,'unc.'), td.thresDesc = 'none'; end
u   = str2double(td.u);
str = td.thresDesc;
