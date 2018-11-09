function DCM = spm_dcm_load(P, save_mem)
% Load a cell array of DCM filenames into a subjects x models cell array
% FORMAT DCM = spm_dcm_load(P)
%
% P        - a DCM's filename or
%          - a GCM's filename (which contains a cell array of DCM files) or
%          - a cell array of DCM filenames or
%          - a character array of DCM filenames
% save_mem - (Optional) if true, only loads priors, posteriors and F for
%            models 2-N
%
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% Peter Zeidman
% $Id: spm_dcm_load.m 7479 2018-11-09 14:17:33Z peter $

if ischar(P)    
    if size(P,1) > 1
        % A character array was provided
        P = cellstr(P);
    else
        % A GCM / DCM filename was provided
        x = load(P);    
        if isfield(x,'GCM')
            P = x.GCM;
        elseif isfield(x,'DCM')
            P = {x.DCM};
        else
            error('Unknown DCM file type');
        end
    end
end

if nargin < 2
    save_mem = false;
end

DCM = cell(size(P));

fprintf('Loading DCMs...');
for s = 1:size(P,1)
    for m = 1:size(P,2)                
        
        if isstruct(P{s,m})
            % DCM already a structure - leave unchanged
            dcm = P{s,m};
        else
            % Load DCM from filename
            try
                model = load(P{s,m});
                dcm   = model.DCM;                
            catch
               fprintf('File: %s\n',P{s,m});
               error('Failed to load model for subject %d model %d', s, m);                
            end
        end
                                
        if save_mem && (m>1)
            % Onlt import basic structures
            DCM{s,m} = struct();
            if isfield(dcm,'a'),  DCM{s,m}.a  = dcm.a; end
            if isfield(dcm,'b'),  DCM{s,m}.b  = dcm.b; end
            if isfield(dcm,'c'),  DCM{s,m}.c  = dcm.c; end
            if isfield(dcm,'d'),  DCM{s,m}.d  = dcm.d; end
            if isfield(dcm,'F'),  DCM{s,m}.F  = dcm.F; end
            if isfield(dcm,'Ep'), DCM{s,m}.Ep = dcm.Ep; end
            if isfield(dcm,'Cp'), DCM{s,m}.Cp = dcm.Cp; end
            if isfield(dcm,'options'), DCM{s,m}.options = dcm.options; end
            if isfield(dcm,'M')
                if isfield(dcm.M,'pE'), DCM{s,m}.M.pE = dcm.M.pE; end
                if isfield(dcm.M,'pC'), DCM{s,m}.M.pC = dcm.M.pC; end
            end
        else
            % Import the whole DCM
            DCM{s,m} = dcm;
        end
               
    end
end
fprintf('Done\n');
