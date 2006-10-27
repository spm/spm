function DCM = spm_dcm_eeg_selectdata(DCM, Spatialmodel, Nmodes)

if Spatialmodel == 1
    % SVD
    
    % concatenate data
    y = [];
    for i = 1:length(DCM.xY.xy)
        y = [y DCM.xY.xy{i}'];
    end

    % maybe include that you use a subset of time points only for data
    % selection -> are there other options for data selection?

    [u s v] = svd(y', 0);
    E = (y*u(:, 1:Nmodes)*s(1:Nmodes, 1:Nmodes)^(-1))';

    % projection of data
    for i = 1:length(DCM.xY.xy)
        DCM.xY.xy{i} = DCM.xY.xy{i}*E';
    end
else
    % no selection, keep all data
    E = eye(size(DCM.xY.xy{1}, 1));
end

DCM.E = E;