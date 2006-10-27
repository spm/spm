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
        P.Tselection = [DCM.xY.Time(1) DCM.xY.Time(end)];
    end
    
    [m, T0] = min(abs(DCM.xY.Time-P.Tselection(1)));
    [m, T1] = min(abs(DCM.xY.Time-P.Tselection(2)));

    % concatenate data
    y = [];
    for i = 1:length(DCM.xY.xy)
        y = [y DCM.xY.xy{i}(T0:T1, :)'];
    end

    % maybe include that you use a subset of time points only for data
    % selection -> are there other options for data selection?

    [u s v] = svd(y', 0);
    E = (y*u(:, 1:P.Nmodes)*s(1:P.Nmodes, 1:P.Nmodes)^(-1))';

    % projection of data
    for i = 1:length(DCM.xY.xy)
        DCM.xY.xy{i} = DCM.xY.xy{i}*E';
    end
    
elseif P.projection == 2
    % dipole leadfields
    D = P.D;
    
    M = DCM.M;
    
    % transformation matrix from MNI-oriented coordinate system
    % to fieldtrip
    iMt = [[0 1 0 20]; [-1 0 0 0]; [0 0 1 10]; [0 0 0 1]];
    iSt = [[0 1 0 0]; [-1 0 0 0]; [0 0 1 0]; [0 0 0 1]];
    P.Lpos = iMt*[P.Lpos; ones(1, size(P.Lpos, 2))];
    P.Lmom = iSt*[P.Lmom; ones(1, size(P.Lmom, 2))];
    P.Lpos = P.Lpos(1:3, :); P.Lmom = P.Lmom(1:3, :);

    for i = 1:M.Nareas
        Lf = fieldtrip_eeg_leadfield4(P.Lpos(:,i), M.dipfit.elc, M.dipfit.vol);
        L(:,i) = Lf*P.Lmom(:,i);
    end

    % normalise
%     L = L./repmat(sqrt(sum(L.^2)), size(L, 1), 1);
     L = L*2*10^4;
    
    E = pinv(L);
    
    % projection of data
    for i = 1:length(DCM.xY.xy)
        DCM.xY.xy{i} = DCM.xY.xy{i}*E';
    end

    % display
    CTF = load(fullfile(spm('dir'), 'EEGtemplates', D.channels.ctf));
    CTF.Cpos = CTF.Cpos(:, D.channels.order(D.channels.eeg));

    x = min(CTF.Cpos(1,:)):0.005:max(CTF.Cpos(1,:));
    y = min(CTF.Cpos(2,:)):0.005:max(CTF.Cpos(2,:));

    [x1, y1] = meshgrid(x, y);
    xp = CTF.Cpos(1,:)';
    yp = CTF.Cpos(2,:)';

    Nchannels = size(CTF.Cpos, 2);

    for i = 1:M.Nareas

        xs = zeros(M.Nareas, 9);
        xs(i, 9) = 1;
        xs = xs(:);
        xs = [0; xs];
    
        ym = NaN*ones(1, Nchannels);
        ym(:, DCM.M.dipfit.chansel) = L(:, i);
        
        figure
        z = griddata(xp, yp, ym, x1, y1);
        surface(x, y, z);
        axis off
        axis square

        shading('interp')
        hold on
        plot3(xp, yp, ym, 'k.');
        title(DCM.Sname{i});

    end

    
else
    % no selection, keep all data
    E = eye(size(DCM.xY.xy{1}, 2));
end

DCM.E = E;