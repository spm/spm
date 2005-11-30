function [c,L,j,res] = spm_eeg_ip_costSd(loc) 

% Cost function for the dipoles fitting routine.
% A "realistic sphere" model is used.
%
% 'loc' represents the location of the n_dip dipoles considered 
% in a 3xn_dip matrix. It is expressed in mm !!!
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Christophe Phillips,
% $Id$

global MODEL V_EEG V_BR OR_OPT FXD_OR

% global variables to avoid passing them each time:
% - MODEL:  contains all info about model (head surfaces, electrodes, sigmas, IFS)
%           'model' was modified in order to include surface vertices coord in mm
%            + electrodes locations in mm
%            + weights for forward problem solution
% - V_EEG:      EEG data to be fitted
% - V_BR:   brain mask to limit the dipoles location.
% - OR_OPT: option for the orientation of the sources
%               1   free orientation.
%               2   fixed orientation, weighted (EEG power) mean of orientations over time.
%               3   fixed orientation, as at the max of EEG power.
%               4   fixed orientation, as defined by the user.

if exist('OR_OPT')==0
    OR_OPT = 1;
end

% m_dist = zeros(MODEL.n_dip,1) ; % mininmum distance between dipoles and brain surface vertices
% br     = zeros(MODEL.n_dip,1) ; % Brain (.i.e brain mask) value at dipoles

m_dist = MODEL.spheres.Rbr-sqrt(sum(loc.^2)) ; % mininmum distance between dipoles and brain surface
if MODEL.n_dip>1
    dist_inter_dip = zeros((MODEL.n_dip-1)*MODEL.n_dip/2,1);
    for ii=1:(MODEL.n_dip-1)
        for jj=(ii+1):MODEL.n_dip
            dist_inter_dip( (ii-1)*(2*MODEL.n_dip-ii)/2 + jj-ii ) = norm(loc(:,ii)-loc(:,jj));
        end
    end
end

L = spm_eeg_inv_Lana(loc,MODEL.spheres.Sc_elXYZ, ...
            MODEL.spheres.Rsc,MODEL.spheres.Rsk,MODEL.spheres.Rbr,MODEL.sigma);
L = L - ones(MODEL.electrodes.nr,1)*mean(L);

% Find "real" location of dipoles (Rloc), to check it's inside the BrainMask!
[Rloc] = spm_eeg_inv_vert2Rsph(2,MODEL.spheres,loc);
Rloc_vx = V_BR.mat\[Rloc ; ones(1,MODEL.n_dip)];
br = spm_sample_vol(V_BR,Rloc_vx(1,:)',Rloc_vx(2,:)',Rloc_vx(3,:)',1);

penal_d = 2 ;
decr_d = 5 ;
penal_out_br = 5;
mult_weight2 = prod(exp(-m_dist/decr_d)*penal_d+1);
    % dipoles to close to the surface are penalized
% mult_weight1 = 1;
mult_weight1 = prod((1-br+eps)*penal_out_br+1);
    % For any dipole outside the brain volume, cost is multiplied by 'penal_out_br'
if MODEL.n_dip>1
    mult_weight3 = prod(2*exp(-(dist_inter_dip-5)/5)+1);
else
    mult_weight3 = 1;
end

% disp([mult_weight1   mult_weight2 mult_weight3 ])

% Sources orientation is fixed according to OR_OPT
%               1   free orientation.
%               2   fixed orientation, weighted (EEG power) mean of orientations over time.
%               3   fixed orientation, as at the max of EEG power.
%               4   fixed orientation, as defined by the user
if OR_OPT==1
    j = pinv(L)*V_EEG;
elseif OR_OPT==2
    j = reshape(pinv(L)*V_EEG,3,MODEL.n_dip*MODEL.Ntb);
    or = j./( ones(3,1) * sqrt(sum(j.^2)) )*MODEL.orW;
    or = or./( ones(3,1) * sqrt(sum(or.^2)) );
    Or = kron(diag(or(1,:)'),[1 0 0]') + kron(diag(or(2,:)'),[0 1 0]') + kron(diag(or(3,:)'),[0 0 1]');
    j = (Or*pinv(L*Or))*V_EEG;   
elseif OR_OPT==3
    jM = reshape(pinv(L)*V_EEG(:,MODEL.Mtb),3,MODEL.n_dip);
    or = jM./( ones(3,1) * sqrt(sum(jM.^2)) );
    Or = kron(diag(or(1,:)'),[1 0 0]') + kron(diag(or(2,:)'),[0 1 0]') + kron(diag(or(3,:)'),[0 0 1]');
    j = (Or*pinv(L*Or))*V_EEG;
elseif OR_OPT==4
    Or = kron(diag(FXD_OR(1,:)'),[1 0 0]') + kron(diag(FXD_OR(2,:)'),[0 1 0]') ...
                    + kron(diag(FXD_OR(3,:)'),[0 0 1]');   
    j = (Or*pinv(L*Or))*V_EEG;
end

res = norm(V_EEG-L*j, 'fro') ; % Fitting of dipole(s) on EEG data
c = res * mult_weight1 * mult_weight2  * mult_weight3;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
