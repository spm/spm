function [L] = spm_erp_L(P,M)
% returns [projected] lead field L as a function of position and momments
% FORMAT [L] = spm_erp_L(P,M)
% P  - model parameters
% M  - model specification
% L  - lead field
%__________________________________________________________________________
% For ECD solutions:
% This routine will autmatically place the lead field in M.L and record the
% paramters in M.Lpos and M.Lmon.  This enables spm_erp_L to compute a new
% lead field only when necessary (i.e., when the parameters change)
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
try
    n = size(P.Lpos,2);
catch
    try
        n = size(P.L,2);
    catch
        n = length(M.pE.A{1});
    end
end

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
        %--------------------------------------------------------------
        LastLpos = P.Lpos;

        % note that in all DCM code, EEG coordinates are in MNI
        % space-orientation. Transform here to POL space and add
        % translation to relate coordinates to centre of sphere
        %--------------------------------------------------------------
        iMt  = M.dipfit.Mmni2polsphere;
        iSt  = iMt; iSt(1:3,4) = 0;

        Lpos = iMt*[P.Lpos; ones(1,n)];
        Lmom = iSt*[P.Lmom; ones(1,n)];
        Lpos = Lpos(1:3,:);
        Lmom = Lmom(1:3,:);

        for i = Id
            Lf = fieldtrip_eeg_leadfield4(Lpos(:,i), M.dipfit.elc, M.dipfit.vol);
            LastL(:,:,i) = Lf*1000;
        end
        for i = 1:n
            L(:,i) = LastL(:,:,i)*Lmom(:,i);
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

        for i = Id
            Lf = fieldtrip_meg_leadfield(P.Lpos(:,i)', M.grad, M.dipfit.vol);
            LastL(:,:,i) = Lf(M.dipfit.Ic,:)*1000;
        end
        for i = 1:n
            L(:,i) = LastL(:,:,i)*P.Lmom(:,i);
        end

    % Imaging solution {specified in M.dipfit.G}
    %----------------------------------------------------------------------
    case{'Imaging'}
        
        % re-compute lead field only if any coeficients have changed
        %------------------------------------------------------------------
        try
            Id = find(any(LastLpos ~= P.L));
        catch
            Id = 1:n;
        end
        for i = Id
            LastL(:,i) = M.dipfit.G{i}*P.L(:,i);
        end
        
        % record new spatial parameters
        %--------------------------------------------------------------
        LastLpos = P.L;
        L        = LastL;

    % LFP or sources - no lead feild required
    %----------------------------------------------------------------------
    case{'LFP'}

        L = sparse(diag(P.L)); return

    otherwise
        warndlg('unknown type of lead field')
end

