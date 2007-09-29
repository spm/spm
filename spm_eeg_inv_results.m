function [D] = spm_eeg_inv_results(D)
% contrast of evoked responses and power for an MEG-EEG model
% FORMAT [D] = spm_eeg_inv_results(D)
% Requires:
%
%     D.inv{i}.contrast.woi   - time (ms) window of interest
%     D.inv{i}.contrast.fboi  - frequency window of interest
%
% this routine will create a contrast for each trial type and will compute
% induced repeonses in terms of power (over trials)
%__________________________________________________________________________

% SPM data structure
%==========================================================================
fprintf('\ncomputing contrast - please wait\n')
try
    model = D.inv{D.val};
catch

    model = D.inv{end};
end

% defaults
%--------------------------------------------------------------------------
try, woi = model.contrast.woi;  catch, woi = [80 120]; end
try, foi = model.contrast.fboi; catch, foi = 0;        end
try, Han = model.contrast.Han;  catch, Han = 1;        end

% inversion parameters
%--------------------------------------------------------------------------
J    = model.inverse.J;                           % Trial average MAP estimate
T    = model.inverse.T;                           % temporal projector
U    = model.inverse.U;                           % spatial  projector
R    = model.inverse.R;                           % referencing matrix
Is   = model.inverse.Is;                          % Indices of ARD vertices
Ic   = model.inverse.Ic;                          % Indices of channels
It   = model.inverse.It;                          % Indices of time bins
pst  = model.inverse.pst;                         % preistimulus tim (ms)
Nd   = model.inverse.Nd;                          % number of mesh dipoles
Nb   = size(T,1);                                 % number of time bins
Nc   = size(U,1);                                 % number of channels

% time-frequency contrast
%==========================================================================

% get [Gaussian] time window
%--------------------------------------------------------------------------
toi  = round(woi*(D.Radc/1000)) + D.events.start + 1;
if toi(1) < It(1) | toi(2) > It(end)
     errordlg('Contrast outside range of sampled window');
     return;
else
     toi = toi - It(1) + 1;
end

if Han
    fwhm = min(diff(toi),8);
    t    = exp(-4*log(2)*([1:Nb] - mean(toi)).^2/(fwhm^2));
    t    = t/sum(t);
else
    t    = toi(1):toi(end);
    t    = ones(size(t))/length(t);
end

% get frequency space and put PST subspace into contrast (W -> T*T'*W)
%--------------------------------------------------------------------------
if foi
    wt = 2*pi*[1:Nb]'/D.Radc;
    W  = [];
    for f = foi(1):foi(end)
        W = [W sin(f*wt) cos(f*wt)];
    end
    W  = diag(t)*W;
    W  = spm_svd(W,1);
else
    W  = t(:);
end
TW     = T'*W;
TTW    = T*TW;

% cycle over trial types
%==========================================================================
trial = model.inverse.trials;
for i = 1:length(J)

    % windowed time average
    %----------------------------------------------------------------------
    if foi

        JW{i} = J{i}*TW(:,1);
        GW{i} = sum((J{i}*TW).^2,2);

    % energy over frequencies
    %----------------------------------------------------------------------
    else

        % get projectors, inversion and data
        %==================================================================
        M     = model.inverse.M;               % MAP projector
        V     = model.inverse.qV;              % temporal correlations
        MUR   = M*U'*R;
        qC    = model.inverse.qC*trace(TTW'*V*TTW);
        
        JW{i} = sparse(0);
        JWWJ  = sparse(0);
        if isfield(D.events,'reject')
            c = find(D.events.code == trial(i) & ~D.events.reject);
        else
            c = find(D.events.code == trial(i));
        end

        % conditional expectation of contrast (J*W) and its energy
        %------------------------------------------------------------------
        Nt    = length(c);
        for j = 1:Nt
            fprintf('\nevaluating trial %i, condition %i',j,i)
            MYW   = MUR*(D.data(Ic,It,c(j))*TTW)/Nt;
            JW{i} = JW{i} + MYW(:,1);
            JWWJ  = JWWJ  + sum(MYW.^2,2);
        end


        % conditional expectation of total energy (source space GW)
        %------------------------------------------------------------------
        GW{i}   = JWWJ + qC/Nt;

        % conditional expectation of induced energy (source space IW)
        % NB: this is zero if Nt = 1
        %--------------------------------------------------------------------------
        % IW   = GW - sum(JW.^2,2) - qC/(Nt(i)*Nt(i));

    end

end



% Save results
%==========================================================================
model.contrast.woi  = woi;
model.contrast.fboi = foi;

model.contrast.W    = W;
model.contrast.JW   = JW;
model.contrast.GW   = GW;

D.inv{D.val}        = model;


% Display
%==========================================================================
spm_eeg_inv_results_display(D);


