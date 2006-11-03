function [D] = spm_eeg_results(D)
% contrast of evoked responses and power for an MEG-EEG model
% FORMAT [D] = spm_eeg_results(D)
% Requires:
%
%     D.contrast.woi   - time (ms) window of interest
%     D.contrast.fboi  - frequency window of interest
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

% inversion parameters
%--------------------------------------------------------------------------
con  = model.inverse.con;                         % condition
M    = model.inverse.MAP;                         % MSP projector
qC   = model.inverse.qC;                          % spatial covariances
V    = model.inverse.qV;                          % temporal correlations
T    = model.inverse.T;                           % PST subsapce
Is   = model.inverse.Is;                          % Indices of ARD vertices
pst  = model.inverse.pst;                         % preistimulus tim (ms)
Nd   = model.inverse.Nd;                          % number of mesh dipoles
Nb   = D.Nsamples;                                % number of time bins
Nc   = size(M,2);                                 % number of channels

% time-frequency contrast
%==========================================================================

% get [Gaussian] time window
%--------------------------------------------------------------------------
toi  = round(woi*(D.Radc/1000)) + D.events.start + 1;
fwhm = diff(toi);
t    = exp(-4*log(2)*([1:Nb] - mean(toi)).^2/(fwhm^2));
t    = t/sum(t);

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
TTW    = T*(T'*W);

% get inversion and data
%==========================================================================

% single-trial analysis (or single ERP)
%--------------------------------------------------------------------------
j     = setdiff(D.channels.eeg, D.channels.Bad);
JW    = sparse(0);
JWWJ  = sparse(0);
c     = find(D.events.code == D.events.types(con));
Nt    = length(c);
TTW   = TTW/Nt;
for i = 1:Nt

    % conditional expectation of contrast (J*W) and its energy
    %-------------------------------------------------------------------
    MYW  = M*squeeze(D.data(j,:,c(i)))*TTW;
    JW   = JW + MYW;
    JWWJ = JWWJ + sum(MYW.^2,2);

end

% conditional expectation of total energy (source space GW)
%--------------------------------------------------------------------------
QC   = qC*trace(TTW'*V*TTW);
GW   = JWWJ + QC*Nt;

% conditional expectation of induced energy (source space IW)
%--------------------------------------------------------------------------
IW   = GW - sum(JW.^2,2) - QC;

% Save results
%==========================================================================
model.contrast.woi    = woi;
model.contrast.fboi   = foi;

model.contrast.W      = W;
model.contrast.JW     = JW;
model.contrast.GW     = GW;
model.contrast.IW     = IW;

D.inv{D.val}         = model;
if str2num(version('-release')) >= 14
    save(fullfile(D.path,D.fname), '-V6','D');
else
    save(fullfile(D.path,D.fname),'D');
end

% Display
%==========================================================================
spm_eeg_results_display(D);


