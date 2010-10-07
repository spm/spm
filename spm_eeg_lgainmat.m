function [L,D] = spm_eeg_lgainmat(D,Is, channels)
% loads or computes if necessary a gain matrix
% FORMAT [L,D] = spm_eeg_lgainmat(D,Is)
% D    - Data structure
% Is   - indices of vertices
%
% L    - Lead-field or gain matrix L(:,Is)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_eeg_lgainmat.m 4082 2010-10-07 16:11:48Z guillaume $


% get gain or lead-field matrix
%--------------------------------------------------------------------------
val = D.val;

forward = D.inv{val}.forward;

for ind = 1:numel(forward)
    modality = forward(ind).modality;
    
    % channels
    %----------------------------------------------------------------------
    chanind = strmatch(modality, D.chantype);
    chanind = setdiff(chanind, D.badchannels);

    if ~isempty(chanind)
        forward(ind).channels = D.chanlabels(chanind);
    else
        error(['No good ' modality ' channels were found.']);
    end
end

if nargin < 3
    channels = [forward(:).channels];
end

try
    fname = D.inv{val}.gainmat;
    G = load(fullfile(D.path, fname)); % Relative path
    
    label = G.label;
    G     = G.G;
    if numel(label) ~= size(G, 1) || ~all(ismember(channels, label))
        error('Gain matrix has an incorrect number of channels');
    end
catch
    label = {};
    for ind = 1:numel(forward)
        % create a new lead-field matrix
        %------------------------------------------------------------------

        % Head Geometry (create tesselation file)
        %------------------------------------------------------------------
        vert = forward(ind).mesh.vert;
        face = forward(ind).mesh.face;

        % normals
        %------------------------------------------------------------------
        norm = spm_mesh_normals(struct('faces',face,'vertices',vert),true);

        vol  = forward(ind).vol;
        
        if ischar(vol)
            vol = ft_read_vol(vol);
        end
        
        sens = D.inv{val}.datareg(ind).sensors;
        modality = forward(ind).modality;     

        % Forward computation
        %------------------------------------------------------------------
        [vol, sens] = ft_prepare_vol_sens(vol, sens, 'channel', forward(ind).channels);
        nvert = size(vert, 1);

        spm('Pointer', 'Watch');drawnow;
        spm_progress_bar('Init', nvert, ['Computing ' modality ' leadfields']); drawnow;
        if nvert > 100, Ibar = floor(linspace(1, nvert,100));
        else Ibar = [1:nvert]; end

        Gxyz = zeros(length(forward(ind).channels), 3*nvert);
        for i = 1:nvert

            Gxyz(:, (3*i- 2):(3*i))  = ft_compute_leadfield(vert(i, :), sens, vol);

            if ismember(i, Ibar)
                spm_progress_bar('Set', i); drawnow;
            end

        end

        spm_progress_bar('Clear');
        spm_progress_bar('Init', nvert, ['Orienting ' modality ' leadfields']); drawnow;

        G{ind} = zeros(size(Gxyz, 1), size(Gxyz, 2)/3);
        for i = 1:nvert

            G{ind}(:, i) = Gxyz(:, (3*i- 2):(3*i))*norm(i, :)';

            if ismember(i, Ibar)
                spm_progress_bar('Set', i); drawnow;
            end

        end

        spm_progress_bar('Clear');

        spm('Pointer', 'Arrow');drawnow;
        
        
        label = [label; forward(ind).channels(:)];
    end
    
    if numel(G)>1
        G = cat(1, G{:});
    else
        G = G{1};
    end
    
    % Save
    %----------------------------------------------------------------------
    D.inv{val}.gainmat = ['SPMgainmatrix_' spm_str_manip(D.fname, 'tr') '_' num2str(val) '.mat'];
    save(fullfile(D.path, D.inv{val}.gainmat), 'G', 'label');
    save(D);
end

[sel1, sel2] = spm_match_str(channels, label);

if length(sel2) ~= numel(channels)
    error('Did not find a match for all the requested channels');
end

% condition the scaling of the lead-field
%--------------------------------------------------------------------------
L   = spm_cond_units(sparse(G(sel2, :)));


% retain selected sources if necessary
%--------------------------------------------------------------------------
if nargin > 1 && ~isempty(Is)
    L = L(:,Is);
end

D.inv{val}.forward = forward;
