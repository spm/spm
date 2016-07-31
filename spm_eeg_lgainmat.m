function [L,D] = spm_eeg_lgainmat(D,Is,channels)
% loads or computes if necessary a gain matrix
% FORMAT [L,D] = spm_eeg_lgainmat(D,Is)
% D    - Data structure
% Is   - indices of vertices
%
% L    - Lead-field or gain matrix L(:,Is)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_eeg_lgainmat.m 6849 2016-07-31 12:34:33Z karl $


% get gain or lead-field matrix
%--------------------------------------------------------------------------
val = D.val;

forward = D.inv{val}.forward;

for ind = 1:numel(forward)
    modality = forward(ind).modality;
    
    % channels
    %----------------------------------------------------------------------
    if isequal(modality, 'MEG')
        chanind = D.indchantype({'MEG', 'MEGPLANAR'}, 'GOOD');
    else
        chanind = D.indchantype(modality, 'GOOD');
    end
    
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
    G     = {};
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
                
        modality = forward(ind).modality;
        
        if isfield(forward, 'siunits') && forward(ind).siunits
            units = D.units(D.indchannel(forward(ind).channels));
            sens  = forward(ind).sensors;
            siunits = isempty(strmatch('unknown', units));           
        else
            siunits = false;
            sens = D.inv{val}.datareg(ind).sensors;
        end            
        
        % Forward computation
        %------------------------------------------------------------------
        [vol, sens] = ft_prepare_vol_sens(vol, sens, 'channel', forward(ind).channels);
        nvert = size(vert, 1);
        
        spm('Pointer', 'Watch');drawnow;
        spm_progress_bar('Init', nvert, ['Computing ' modality ' leadfields']); drawnow;
        if nvert > 100, Ibar = floor(linspace(1, nvert,100));
        else Ibar = [1:nvert]; end
        
        if ~isequal(ft_voltype(vol), 'interpolate')
            Gxyz = zeros(length(forward(ind).channels), 3*nvert);
            for i = 1:nvert
                if siunits
                    Gxyz(:, (3*i - 2):(3*i))  = ft_compute_leadfield(vert(i, :), sens, vol,...
                        'dipoleunit', 'nA*m', 'chanunit', units);
                else
                    Gxyz(:, (3*i - 2):(3*i))  = ft_compute_leadfield(vert(i, :), sens, vol);
                end
                
                if ismember(i, Ibar)
                    spm_progress_bar('Set', i); drawnow;
                end
            end
        else
            if siunits
                Gxyz = ft_compute_leadfield(vert, sens, vol, 'dipoleunit', 'nA*m', 'chanunit', units);
            else
                Gxyz = ft_compute_leadfield(vert, sens, vol);
            end
        end
        
        spm_progress_bar('Clear');
        spm_progress_bar('Init', nvert, ['Orienting ' modality ' leadfields']); drawnow;
        
        G{ind} = zeros(size(Gxyz, 1), size(Gxyz, 2)/3);
        for i = 1:nvert
            G{ind}(:, i) = Gxyz(:, (3*i- 2):(3*i))*norm(i, :)';
            if ismember(i,Ibar)
                spm_progress_bar('Set', i); drawnow;
            end
            
        end
        
        % condition the scaling of the lead-field
        %--------------------------------------------------------------------------        
        [Gs, scale] = spm_cond_units(G{ind});
     
        if siunits && abs(log10(scale))>2
            warning(['Scaling expected to be 1 for SI units, actual scaling ' num2str(scale)]);
            G{ind} = Gs;
        else
            scale = 1;
        end
        
        spm_progress_bar('Clear');
        
        spm('Pointer', 'Arrow');drawnow;
        
        label = [label; forward(ind).channels(:)];
        
        forward(ind).scale = scale;
    end
    
    if numel(G) > 1
        G = cat(1, G{:});
    else
        G = G{1};
    end
    
    % Save
    %----------------------------------------------------------------------
    D.inv{val}.gainmat = ['SPMgainmatrix_' spm_file(D.fname, 'basename') '_' num2str(val) '.mat'];
    save(fullfile(D.path, D.inv{val}.gainmat), 'G', 'label');
    save(D);
end

[sel1, sel2] = spm_match_str(channels, label);

if length(sel2) ~= numel(channels)
    error('Did not find a match for all the requested channels');
end

L   = sparse(G(sel2, :));

% retain selected sources if necessary
%--------------------------------------------------------------------------
if nargin > 1 && ~isempty(Is)
    L = L(:,Is);
end

D.inv{val}.forward = forward;
