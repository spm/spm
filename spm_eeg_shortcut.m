function [Ishortcut, SPM] = spm_eeg_shortcut(SPM)
% function SPM = spm_eeg_shortcut(SPM)
% Internally used function that shortcuts computation of SVD of simple EEG
% design matrices
% Output:
%    Ishortcut: 0/1 indicator whether shortcut could be taken
%    SPM: SPM struct with updated elements if shortcut was taken
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_shortcut.m 213 2005-08-22 12:43:29Z stefan $

Ishortcut = 0;

X1 = speye(size(SPM.xX.X, 1));

if all(size(SPM.xX.X) == size(X1)) 
    if SPM.xX.X == X1
        % big identity matrix
        Ishortcut = 1;
        
        SPM.Vbeta = SPM.xY.VY;
        SPM.xVol.M = SPM.xY.VY(1).mat;
        SPM.xVol.DIM = SPM.xY.VY(1).dim(1:3)';
        
        SPM.Vbeta = SPM.xY.VY;
        SPM.xVol.M = SPM.xY.VY(1).mat;
        SPM.xVol.DIM = SPM.xY.VY(1).dim(1:3)';
        
        % instead of: 
        % SPM.xX.xKXs = spm_sp('Set', spm_filter(1, SPM.xX.X));		% KWX
        n = size(SPM.xX.X, 1);
        xKXs.X = SPM.xX.X;
        xKXs.ds = ones(n, 1);
        xKXs.tol =  n*max(abs(xKXs.ds))*eps;
        xKXs.u = speye(n);
        xKXs.v = speye(n);
        xKXs.rk = n;
        xKXs.oP = [];
        xKXs.oPp = [];
        xKXs.ups = [];
        xKXs.sus = [];
        
        SPM.xX.xKXs = xKXs;
        
        SPM.xCon = [];
    end
end