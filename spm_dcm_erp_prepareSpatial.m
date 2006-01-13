function DCM = spm_dcm_erp_prepareSpatial(DCM)
% prepares structures for ECD forward model (both EEG and MEG)

% Stefan Kiebel
% $Id: spm_dcm_erp_prepareSpatial.m 404 2006-01-13 18:42:21Z stefan $

if DCM.options.Spatial_type == 1
    % EEG
    
    sensorfile = DCM.M.dipfit.sensorfile;
    
    % Polyhmus file
    if strcmpi('pol', spm_str_manip(sensorfile, 'e'))
        try
            [x, y, z] = textread(sensorfile, '%f %f %f', 'headerlines', 6);
            x = x(5:end-1); y = y(5:end-1); z = z(5:end-1);
        catch
            errordlg('Could not read sensor location file');
        end
    else
        try
            tmp = load(sensorfile, '-mat');
            xyz = fieldnames(tmp);
            eval(['xyz = tmp.' xyz{1}]);
           
            x = xyz(:, 1); y = xyz(:, 2); z = xyz(:, 3);
        catch
            errordlg('Could not read sensor location file');
        end
    end

    dipfit = DCM.M.dipfit;

    dipfit.chansel = DCM.M.Ichannels; % remove bad channels
    dipfit.elc = [x(dipfit.chansel) y(dipfit.chansel) z(dipfit.chansel)];

    dipfit.vol.r = [71 72 79 85];
    dipfit.vol.c = [0.3300 1 0.0042 0.3300];
    dipfit.vol.o = [0 0 0];

    melc = mean(dipfit.elc);
    dipfit.elc = dipfit.elc - repmat(melc, size(dipfit.elc, 1), 1);

    % projecting to outer sphere
    dist = sqrt(sum(dipfit.elc.^2,2));
    dipfit.elc = dipfit.vol.r(4) * dipfit.elc ./ [dist dist dist];

    DCM.M.dipfit = dipfit;
elseif DCM.options.Spatial_type == 2
    % MEG
    
    % not much to do, because the sensors location/orientations were
    % already read at the time of conversion. 
    
    DCM.M.dipfit.vol.r = [85];
    DCM.M.dipfit.vol.o = [0 0 0];
end
