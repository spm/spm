function [L] = spm_sep_L(P,M)
% returns [projected] lead field L as a function of position and moments
% FORMAT [L] = spm_sep_L(P,M)
% P  - model parameters
% M  - model specification
% L  - lead field
%__________________________________________________________________________
% For ECD solutions:
% This routine will automatically place the lead field in M.L and record the
% parameters in M.Lpos and M.Lmon.  This enables spm_sep_L to compute a new
% lead field only when necessary (i.e., when the parameters change)
%
% This routine is the same as spm_erp_L, with the expectation that two sets
% of lead fields are produced; one for the spiny stellate cells and the
% other for pyramidal cells
%
% see; Kiebel et al. (2006) NeuroImage
%__________________________________________________________________________
% %W% Karl Friston %E%

% Create a persient variable that rembers the last locations
%--------------------------------------------------------------------------
persistent LastLpos LastL

% output
%==========================================================================

% number of sources - n
%--------------------------------------------------------------------------
n = length(M.pE.A{1});

switch M.dipfit.type

    % parameterised lead field ECD/EEG
    %----------------------------------------------------------------------
    case{'ECD (EEG)'}

        % re-compute lead field only if any dipoles changed
        %------------------------------------------------------------------
        try
            Id = find(any(LastLpos ~= P.Lpos));
        catch
            Id = 1:n;
        end

        % record new spatial parameters
        %------------------------------------------------------------------
        LastLpos = P.Lpos;

        % note that in all DCM code, EEG coordinates are in MNI
        % space-orientation. Transform here to POL space and add
        % translation to relate coordinates to centre of sphere
        %------------------------------------------------------------------
        iMt  = M.dipfit.Mmni2polsphere;
        iSt  = iMt; iSt(1:3,4) = 0;

        Lpos = iMt*[P.Lpos; ones(1,n)];
        Lmom = iSt*[P.Lmom; ones(1,n + n)];
        Lpos = Lpos(1:3,:);
        Lmom = Lmom(1:3,:);

        % new lead-feilds
        %------------------------------------------------------------------
        for i = Id
            Lf = fieldtrip_eeg_leadfield4(Lpos(:,i), M.dipfit.elc, M.dipfit.vol);
            LastL(:,:,i) = Lf*20000;
        end
        
        % new moments (pyramidal)
        %------------------------------------------------------------------
        for i = 1:n
            L(:,i) = LastL(:,:,i)*Lmom(:,i);
        end
        
        % new moments (stellate)
        %------------------------------------------------------------------
        for i = 1:n
            L(:,i + n) = LastL(:,:,i)*Lmom(:,i + n);
        end


    % parameterised lead field ECD/MEG
    %----------------------------------------------------------------------
    case{'ECD (MEG)'}

        % re-compute lead field only if any dipoles changed
        %------------------------------------------------------------------
        try
            Id = find(any(LastLpos ~= P.Lpos));
        catch
            Id = 1:n;
        end

        % record new locations and compute new lead field components
        %------------------------------------------------------------------
        LastLpos = P.Lpos;

        % new lead-feilds
        %------------------------------------------------------------------
        for i = Id
            Lf = fieldtrip_meg_leadfield(P.Lpos(:,i)', M.grad, M.dipfit.vol);
            LastL(:,:,i) = Lf(M.dipfit.Ic,:)*1e12;
        end
        
        % new moments (pyramidal)
        %------------------------------------------------------------------
        for i = 1:n
            L(:,i) = LastL(:,:,i)*Lmom(:,i);
        end
        
        % new moments (stellate)
        %------------------------------------------------------------------
        for i = 1:n
            L(:,i + n) = LastL(:,:,i)*Lmom(:,i + n);
        end

    % Imaging solution {specified in M.dipfit.G}
    %----------------------------------------------------------------------
    case{'Imaging'}
        
        % new coeficients (stellate)
        %------------------------------------------------------------------
        for i = 1:n
            L(:,i) = M.dipfit.G{i}*P.L(:,i);
        end
        
        % new coeficients (pyramidal)
        %------------------------------------------------------------------
        for i = 1:n
            L(:,i + n) = M.dipfit.G{i}*P.L(:,i + n);
        end
        

    % LFP or sources - no lead feild required
    %----------------------------------------------------------------------
    case{'LFP'}

        L = sparse(diag(P.L)); return

    otherwise
        warndlg('unknown type of lead field')
end

