function [L] = spm_erp_L(P)
% returns lead field L as a function of position and momments
% FORMAT [L] = spm_erp_L(P)
% P  - model parameters
% L  - lead field
%__________________________________________________________________________
% This routine will autmatically place the lead field in M.L and record the
% paramters in M.Lpos and M.Lmon.  This enables spm_erp_L to compute a new 
% lead field only when necessary (i.e., when the parameters change)
%
% see; Kiebel et al.
%__________________________________________________________________________
% %W% Karl Friston %E%

% output
%==========================================================================
global M

switch M.Spatial_type
        
    % parameterised lead field ECD/EEG
    %----------------------------------------------------------------------
    case{1}

        % re-compute lead field only if any dipoles changed
        %------------------------------------------------------------------
        try
           Id = find(any([M.Lpos; M.Lmom] ~= [P.Lpos; P.Lmom]));
        catch
           Id = 1:size(P.Lpos,2);
        end
        
        % record new spatial parameters
        %--------------------------------------------------------------
        M.Lpos = P.Lpos;
        M.Lmom = P.Lmom;

        if length(Id)

            % note that in all DCM code, EEG coordinates are in MNI
            % space-orientation. Transform here to fieldtrip.
            % transformation matrix from MNI-oriented coordinate system
            % to fieldtrip
            %--------------------------------------------------------------
            iMt = [0 1 0 20;
                  -1 0 0 0;
                   0 0 1 10;
                   0 0 0 1];
            iSt = [0 1 0 0;
                  -1 0 0 0;
                   0 0 1 0;
                   0 0 0 1];
            
            P.Lpos = iMt*[P.Lpos; ones(1, size(P.Lpos, 2))];
            P.Lmom = iSt*[P.Lmom; ones(1, size(P.Lmom, 2))];
            P.Lpos = P.Lpos(1:3,:); P.Lmom = P.Lmom(1:3,:);

            for i = Id
                Lf = fieldtrip_eeg_leadfield4(P.Lpos(:,i), M.dipfit.elc, M.dipfit.vol);
                M.L(:,i) = Lf*P.Lmom(:,i)*20000;
            end
            
        end

    % parameterised lead field ECD/MEG
    %----------------------------------------------------------------------
    case{2}

        % re-compute lead field only if any dipoles changed
        %------------------------------------------------------------------
        try
           Id = find(any([M.Lpos; M.Lmom] ~= [P.Lpos; P.Lmom]));
        catch
           Id = 1:size(P.Lpos,2);
        end
        
        % record new spatial parameters
        %------------------------------------------------------------------
        M.Lpos = P.Lpos;
        M.Lmom = P.Lmom;
            
        if length(Id)
            for i = Id
                Lf = fieldtrip_meg_leadfield(P.Lpos(:,i)', M.grad, M.dipfit.vol);
                M.L(:,i) = Lf(M.Ichannels,:)*P.Lmom(:,i)*(10^20);
            end


        end
        
    % fixed lead field {specified in M.L}
    %----------------------------------------------------------------------
    case{3}

        M.L = M.L;
        
    otherwise
        warndlg('unknown type of lead field')
end

% lead feild
%--------------------------------------------------------------------------
L   = M.L;

