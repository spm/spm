function DCM = spm_dcm_eeg_selectdata(DCM, P)
% function that performs data selection on data in DCM using parameters P
% FORMAT DCM = spm_dcm_eeg_selectdata(DCM, P)
% 
% DCM		    - DCM structure
% (mandatory) fields of P:
% projection    - type of data selection: 0: none, 1: singular value
%                 decomposition
% Nmodes        - number of modes required
% (optional) fields of P:
% Tselection    - start and end in peri-stimulus time. Only data between
%                 these time points is used to compute projection. Default
%                 is all data.
%
% Output:
% DCM			- DCM structure
%_______________________________________________________________________
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel

if P.projection == 1
    % SVD
    
    if ~isfield(P, 'Tselection')
        P.Tselection = [DCM.Y.Time(1) DCM.Y.Time(end)];
    end
    
    [m, T0] = min(abs(DCM.Y.Time-P.Tselection(1)));
    [m, T1] = min(abs(DCM.Y.Time-P.Tselection(2)));

    % concatenate data
    y = [];
    for i = 1:length(DCM.Y.xy)
        y = [y DCM.Y.xy{i}(T0:T1, :)'];
    end

    % maybe include that you use a subset of time points only for data
    % selection -> are there other options for data selection?

    [u s v] = svd(y', 0);
    E = (y*u(:, 1:P.Nmodes)*s(1:P.Nmodes, 1:P.Nmodes)^(-1))';

    % projection of data
    for i = 1:length(DCM.Y.xy)
        DCM.Y.xy{i} = DCM.Y.xy{i}*E';
    end
else
    % no selection, keep all data
    E = eye(size(DCM.Y.xy{1}, 1));
end

DCM.E = E;