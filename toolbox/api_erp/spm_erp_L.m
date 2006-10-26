function [L] = spm_erp_L(P)
% returns [projected] lead field L as a function of position and momments
% FORMAT [L] = spm_erp_L(P)
% P  - model parameters
% L  - lead field
%__________________________________________________________________________
% This routine will autmatically place the lead field in M.L and record the
% paramters in M.Lpos and M.Lmon.  This enables spm_erp_L to compute a new
% lead field only when necessary (i.e., when the parameters change)
%
% see; Kiebel et al. (2006) NeuroImage
%__________________________________________________________________________
% %W% Karl Friston %E%

% output
%==========================================================================
global M
n   = size(P.Lpos,2);                                   % number of sources

switch M.Spatial_type

    % parameterised lead field ECD/EEG
    %----------------------------------------------------------------------
    case{1}

        % re-compute lead field only if any dipoles changed
        %------------------------------------------------------------------
        try
            Id = find(any(M.Lpos ~= P.Lpos));
        catch
            Id = 1:n;
        end

        % record new spatial parameters
        %--------------------------------------------------------------
        M.Lpos = P.Lpos;

        % note that in all DCM code, EEG coordinates are in MNI
        % space-orientation. Transform here to POL space and add
        % translation to relate coordinates to centre of sphere
        %--------------------------------------------------------------
        iMt  = M.dipfit.Mmni2polsphere;
        iSt  = iMt; iSt(1:3,4) = 0;

        Lpos = iMt*[P.Lpos; ones(1, size(P.Lpos, 2))];
        Lmom = iSt*[P.Lmom; ones(1, size(P.Lmom, 2))];
        Lpos = Lpos(1:3,:);
        Lmom = Lmom(1:3,:);

        for i = Id
            Lf = fieldtrip_eeg_leadfield4(Lpos(:,i), M.dipfit.elc, M.dipfit.vol);
            M.L(:,:,i) = Lf*20000;
        end
        for i = 1:n
            L(:,i) = M.L(:,:,i)*Lmom(:,i);
        end


    % parameterised lead field ECD/MEG
    %----------------------------------------------------------------------
    case{2}

        % re-compute lead field only if any dipoles changed
        %------------------------------------------------------------------
        try
            Id = find(any(M.Lpos ~= P.Lpos));
        catch
            Id = 1:n;
        end

        % record new locations and compute new lead field components
        %------------------------------------------------------------------
        M.Lpos = P.Lpos;

        for i = Id
            Lf = fieldtrip_meg_leadfield(P.Lpos(:,i)', M.grad, M.dipfit.vol);
            M.L(:,:,i) = Lf(M.dipfit.chansel,:)*(10^20);
        end
        for i = 1:n
            L(:,i) = M.L(:,:,i)*P.Lmom(:,i);
        end

    % fixed lead field {specified in M.L}
    %----------------------------------------------------------------------
    case{3}

        L = M.L;

    % LFP or sources - no lead feild required
    %----------------------------------------------------------------------
    case{4}

        L = speye(n,n); return

    otherwise
        warndlg('unknown type of lead field')
end

