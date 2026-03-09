function spm_mac_permissions()
% SPM_MAC_PERMISSIONS Helper to fix macOS quarantine issues (Gatekeeper)
%
% This function detects if SPM's MEX files are tainted with the macOS
% "com.apple.quarantine" attribute (common when downloading ZIPs from the
% web). If detected, it attempts to self-repair by removing the attribute.
%
% This function is called automatically at startup by spm.m.
%__________________________________________________________________________
% Copyright (C) 2024 Wellcome Centre for Human Neuroimaging

if ~ismac || isdeployed
    return;
end

% Check a representative MEX file
spmDir = fileparts(which('spm.m'));
mexFiles = dir(fullfile(spmDir, '*.mexmac*'));

if ~isempty(mexFiles)
    % Check if the first found MEX file has the quarantine attribute
    checkCmd = ['xattr -p com.apple.quarantine "' fullfile(mexFiles(1).folder, mexFiles(1).name) '"'];
    [status, result] = system(checkCmd);
    
    % status==0 means attribute exists
    if status == 0 && ~isempty(strtrim(result))
        fprintf('SPM: Detected macOS Quarantine on MEX files. Attempting self-repair...\n');
        
        % Recursively clear attributes from the SPM directory
        [fixStatus, ~] = system(['xattr -cr "' spmDir '"']);
        
        if fixStatus == 0
            fprintf('SPM: Self-repair successful. Quarantine removed.\n');
        else
            fprintf(2, '\n========================================================================\n');
            fprintf(2, 'SPM: Critical Permission Error\n');
            fprintf(2, 'macOS is blocking SPM binaries because they are quarantined.\n');
            fprintf(2, 'Self-repair failed (likely due to permissions).\n');
            fprintf(2, 'Please run the following command in your Terminal to fix this:\n\n');
            fprintf(2, '     sudo xattr -cr "%s"\n\n', spmDir);
            fprintf(2, 'Then restart MATLAB.\n');
            fprintf(2, '========================================================================\n\n');
        end
    end
end

end
