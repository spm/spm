function G = eeg_sph(L,Channel,Param,Order,Verbose,varargin);
%EEG_SPH - Calculate the electric potential , spherical head, arbitrary orientation
% function G = eeg_sph(L,Channel,Param,Order,Verbose,varargin);
% function G = eeg_sph(L,Channel,Param,Order);
% L is 3 x nL, each column a source location
% Channel is the channel structure, same for Param
% Order is 
%  -1 current dipole  
%   0 focal(magnetic) dipole % NOT SUPPORTED
%   1 1st order multipole    % NOT SUPPORTED
% Param is 
%  .EEGType is one of {'EEG_SINGLE', 'EEG_BERG', 'EEG_3SHELL'};
%  .Berg is set in Param(1) as
%    .mu
%    .lam
%  .Radii vector of radii, inside to outside
%  .Conductivity vector of sigmas inside to outside
%  .Center the sphere center
%
% Verbose : toggle Verbose mode
%
% See also BERG

%<autobegin> ---------------------- 27-Jun-2005 10:44:15 -----------------------
% ------ Automatically Generated Comments Block Using AUTO_COMMENTS_PRE7 -------
%
% CATEGORY: Forward Modeling
%
% Alphabetical list of external functions (non-Matlab):
%   toolbox\dlegpoly.m
%   toolbox\dotprod.m
%   toolbox\good_channel.m
%   toolbox\rownorm.m
%
% Subfunctions in this file, in order of occurrence in file:
%   G = gainp_sph6x(Rq,Re,R,sigma,nmax,method,mu_berg_in,lam_berg_in)
%
% At Check-in: $Author: Mosher $  $Revision: 24 $  $Date: 6/27/05 8:59a $
%
% This software is part of BrainStorm Toolbox Version 27-June-2005  
% 
% Principal Investigators and Developers:
% ** Richard M. Leahy, PhD, Signal & Image Processing Institute,
%    University of Southern California, Los Angeles, CA
% ** John C. Mosher, PhD, Biophysics Group,
%    Los Alamos National Laboratory, Los Alamos, NM
% ** Sylvain Baillet, PhD, Cognitive Neuroscience & Brain Imaging Laboratory,
%    CNRS, Hopital de la Salpetriere, Paris, France
% 
% See BrainStorm website at http://neuroimage.usc.edu for further information.
% 
% Copyright (c) 2005 BrainStorm by the University of Southern California
% This software distributed  under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html .
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%<autoend> ------------------------ 27-Jun-2005 10:44:15 -----------------------


% /---Script Authors--------------------------------------\
% |                                                       |
% | *** John Ermer, Ph.D.                                 |
% |    Signal % Image Processing Institute                |
% |    University of Southern California                  |
% |    Los Angeles, CA, USA                               |
% |                                                       |
% | *** John C. Mosher, Ph.D.                             |
% |  Biophysics Group                                    |
% |                                                       |
% | *** Sylvain Baillet Ph.D.                             |
% | Cognitive Neuroscience & Brain Imaging Laboratory     |
% | CNRS UPR640 - LENA                                    | 
% | Hopital de la Salpetriere, Paris, France              |
% | sylvain.baillet@chups.jussieu.fr                      |
% |                                                       |
% \-------------------------------------------------------/
%  
% Date of creation: October, 25 1999
%
% Script History -----------------------------------------------------------------------------------------------------------
% SB  19-Nov-2002 : Edited Header
%                  Updated management of EEG reference
% JCM 20-Nov-2002 : Fixed headers to have only one autocomments block
% SB  09-Mar-2004 : Added verbose mode
% --------------------------------------------------------------------------------------------------------------------------

if nargin < 5
    Verbose = 1; % Default
end

%-----------------------------

nmax = 80;

%-----------------------------

% EEG Channels
EEGndx = good_channel(Channel,[],'EEG');

% Reference Channel
REFndx = good_channel(Channel,[],'EEG REF');

% EEG Channel locations

% Add EEG reference at the end of Channel and Param structures
Re = [Channel(EEGndx).Loc,Channel(REFndx).Loc]'; % Electrode location array
if length(Param) ~= length([EEGndx,REFndx])
    Param(REFndx) = Param(EEGndx(1));
end

Param = Param([EEGndx,REFndx]);

Rq = L';
center = [Param.Center]';
Re = Re - center;

Rq = Rq - repmat(center(1,:),size(Rq,1),1); % Back to origin [0 0 0] for the sensors and the sources

clear tmp 

switch(Param(1).EEGType)
case 'EEG_SINGLE'
    R = Param(1).Radii(end);
    sigma = Param(1).Conductivity(end);
otherwise
    R = Param(1).Radii;
    sigma = Param(1).Conductivity;
end

switch(Param(1).EEGType)
    
case 'EEG_BERG'
    method = 2; 
    mu_berg_in = Param(1).Berg.mu;
    lam_berg_in  = Param(1).Berg.lam;
case 'EEG_3SHELL'   
    method = 1; 
    mu_berg_in = [];
    lam_berg_in = [];
otherwise
    method = 1; 
    mu_berg_in = [];
    lam_berg_in  = [];
end

if 1 % Projection of the EEG sensors on the sphere
    [theta phi Re_sph] = cart2sph(Re(:,1),Re(:,2),Re(:,3));
    Re_sph = R(end)*ones(size(Re_sph));
    [Re(:,1) Re(:,2) Re(:,3)] = sph2cart(theta,phi,Re_sph);
    
end

Gtmp = gainp_sph6x(Rq,Re,R,sigma,nmax,method,mu_berg_in',lam_berg_in');

if isempty(REFndx) % Average Reference
    
    G = NaN * zeros(length(EEGndx),size(Gtmp,2));
    %G(EEGndx,:) = Gtmp - repmat(mean(Gtmp),size(Gtmp,1),1);   
    G = Gtmp - repmat(mean(Gtmp),size(Gtmp,1),1);   
    clear Gtmp
   
else % Specific electrode as reference
    % Remove its lead field from all others
    % Note: reference leadfield is at the end of original G matrix (from gain_sph) 
    G = NaN * zeros(length(EEGndx),size(Gtmp,2));
    %G(EEGndx,:) = Gtmp(1:end-1,:) - repmat(Gtmp(end,:),size(Gtmp,1)-1,1);   
    %G = G(EEGndx,:);
    G = Gtmp(1:end-1,:) - repmat(Gtmp(end,:),size(Gtmp,1)-1,1);   
    
    %     % Reference lead field is zero
    %     G(REFndx,:) = zeros(1,size(G,2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subfunctions -------------------------------------------------------------------------------
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G = gainp_sph6x(Rq,Re,R,sigma,nmax,method,mu_berg_in,lam_berg_in)

%GAINP_SPH6X EEG Multilayer Spherical Forward Model
% function G = gainp_sph6x(Rq,Re,R,sigma,nmax,method,mu_berg_in,lam_berg_in)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EEG MULTILAYER SPHERICAL FORWARD MODEL (gainp_sph6x.m)
%
% This function computes the voltage potential forward gain matrix for an array of 
% EEG electrodes on the outermost layer of a single/multilayer conductive sphere. 
% Each region of the multilayer sphere is assumed to be concentric with 
% isontropic conductivity.  EEG sensors are assumed to be located on the surface
% of the outermost sphere. 
%
% Calculation of the electric potiential is performed using either (user-specified)
% of the following methods (Ref: Z. Zhang "A fast method to compute surface 
% potentials generated by dipoles within multilayer anisotropic spheres" 
% (Phys. Med. Biol. 40, pp335-349,1995)    
%
%   1) Closed Form Solution (Single Shell Case Only). See formulas (1H,1H')
%
%   2) Series Expansion using Legendre Polynomials. See formulas (1I,2I,3I and 4I)
%
%   3) Series Approximiation of a Multilayer Sphere as three dipoles in a 
%      single shell using "Berg/Sherg" parameter approximation.
%      See formulas (1i',5i" and 6i)
%
% Dipole generator(s) are assumed to be interior to the innermost "core" layer. For those 
% dipoles external to the sphere, the dipole "image" is computed and used determine the 
% gain function. The exception to this is the Legendre Method where all dipoles MUST be 
% interior to the innermost "core" layer.
%
% INPUTS (Required):
%       Rq   : dipole location(in meters)                                 P x 3
%       Re   : EEG sensors(in meters) on the scalp                        M x 3
%       R    : radii(in meters) of sphere from 
%              INNERMOST to OUTERMOST                                     NL x 1
%       sigma: conductivity from INNERMOST to OUTERMOST                   NL x 1
%
% INPUTS (Optional):
%       nmax : # of terms used in Truncated Legendre Series               scalar
%              If not specified, a default value based on outermost
%              dipole magnitude is computed. (Note: This parameter
%              is ignored when Berg Approximation is commanded)
%     method : Method used for computing forward potential    
%              1=Legendre Series Approx; 2=Berg Parameter Approx
%              (Note: Default and all other values invoke Legendre 
%              Series Approx. Exception is single-shell case where
%              closed form solution is always used)                       scalar 
%  mu_berg_in: User specified initial value for Berg eccentricity
%              factors (Required if Berg Method is commanded)             3 x 1
% lam_berg_in: User specified initial value for Berg magnitude
%              factors (Required if Berg Method is commanded)             3 x 1
%  
%              WHERE: M=# of sensors; P=# of dipoles; NL = # of sphere layers
%
% OUTPUTS:
%       G    : EEG forward model gain matrix                              M x (3*P)
%
% External Functions and Files:
%    dlegpoly.m; rownorm.m; dotprod.m: USC/LANL MEG/EEG Toolbox
%    zhang_fit.m: External Function used to fit Berg Parameters (Zhang Eq# 5i")
%    
%     - John Ermer 6/3/99 
%        - 8/9/99: Modified to use EM image for dipoles external to brain
%                  (Applies to single-shell and Berg Methods only) (John Ermer)
%        - 10/31/99: Optimized Processing Associated with Dipoles falling outside sphere
%                  (Applies to single-shell and Berg Methods only) (John Ermer)
%        - 01/26/00: Corrected minor dimension error which caused program to fault when 
%                  external dipoles and multiple sensors were present. (John Ermer)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% THIS PART CHECKS INPUT PARAMETERS FOR THEIR DIMENSION AND VALIDITY %%%
%
NL = length(R);                          % # of concentric sphere layers
P = size(Rq,1);
M = size(Re,1);
% 
if R(1)~= min(R)
    error('Head radii must be specified from innermost to outmost layer!!! ')
end
%
if size(Rq,2) ~= 3
    error('Dipole location must have three columns!!!')
end
%
if nargin < 6,     % Check # of input terms to see if method is specified
    method = 1;                        % Default Method = Legendre Series Expansion
else
    if (method>2)|(method<0)
        method = 0;                     % Default Method = Legendre Series Expansion
    end
end
%
%%% This part pre-initializes parameters used in future calculations %%%
%
G = zeros(M,3*P);   % Pre-Allocate Gain Matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This part computes the potential for a dipole contained within a single-layer
%%% homogeneous sphere using closed form formula 
%%% The EEG single-shell solution uses the vector form (which avoids computationally
%%% expensive intrinsic functions) described by Mosher et al ("EEG and MEG: Forward 
%%% Solutions for inverse problems" IEEE BME Trans March 1999)  
%
if NL == 1  % Single Shell Case (Closed Form Solution)
    %
    Re_mag = repmat(R(NL),P,M);           %(PxM)
    Re_mag_sq = repmat(R(NL)*R(NL),P,M);  %(PxM)
    Rq_mag = rownorm(Rq);                 %(Px1)
    %
    Rq1 = Rq;                             %(Px1)
    Rq1_mag = Rq_mag;                     %(Px1)
    Rq1_mag_sq = Rq_mag.*Rq_mag;          %(Px1)
    Re_dot_Rq1 = Rq1*Re';                 %(PxM)
    %
    const = 4.0*pi*sigma(NL);
    term = 1./(const*Rq1_mag_sq);         %(Px1)
    %
    %%% This part checks for the presence of Berg dipoles which are external to
    %%% the sphere. For those dipoles external to the sphere, the dipole parameters
    %%% are replaced with the electrical image (internal to sphere) of the dipole
    %
    nx = find(Rq1_mag > R(NL));
    %
    if nx>0
        Rq1_temp = Rq1(nx,:);
        Rq1(nx,:) = R(NL)*R(NL)*Rq1_temp./repmat((rownorm(Rq1_temp).*rownorm(Rq1_temp)),1,3);
        Rq1_mag(nx,1) = rownorm(Rq1(nx,:));
        Rq1_mag_sq(nx,1) = Rq1_mag(nx,1).*Rq1_mag(nx,1);
        Re_dot_Rq1(nx,:) = R(NL)*R(NL)*Re_dot_Rq1(nx,:)./repmat((rownorm(Rq1_temp).*rownorm(Rq1_temp)),1,M);
        term(nx,:) = (R(NL)./rownorm(Rq1_temp)).*term(nx,:); 
        %  was : term(nx,:) = (R(NL)/rownorm(Rq1_temp))'.*term(nx,:); 
        
    end
    %
    %%% Calculation of Forward Gain Matrix Contribution due to K-th Berg Dipole
    %
    Rq1_mag = repmat(Rq1_mag,1,M);                          %(PxM)
    Rq1_mag_sq = repmat(Rq1_mag_sq,1,M);                    %(PxM)
    term = repmat(term,1,M);                                %(PxM)
    %
    d_mag = reshape( rownorm(reshape(repmat(Re,1,P)',3,P*M)' ...
        -repmat(Rq1,M,1)) ,P,M);           %(PxM)
    d_mag_cub = d_mag.*d_mag.*d_mag;                        %(PxM)
    F_scalar = d_mag.*(Re_mag.*d_mag+Re_mag_sq-Re_dot_Rq1); %(PxM)
    c1 = term.*(2*( (Re_dot_Rq1-Rq1_mag_sq)./d_mag_cub) ...
        + 1./d_mag - 1./Re_mag);                 %(PxM)
    c2 = term.*((2./d_mag_cub) + (d_mag+Re_mag)./(Re_mag.*F_scalar));
    %
    G = G + reshape(repmat((c1 - c2.*Re_dot_Rq1)',3,1),M,3*P) ...
        .*repmat(reshape(Rq1',1,3*P),M,1) ...
        +  reshape(repmat((c2.*Rq1_mag_sq)',3,1),M,3*P) ...
        .*repmat(Re,1,P);
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %%% This part computes the potential for a dipole contained within a multi-layer
    %%% isontropic sphere using Legendre Polynomial Expansion (Zhang Eqs 1H, 1H')  
    %%% (Code based on gainp_sph.m by CCH,  Aug/20/1995)
    %
elseif (NL>1)&(method==1)  % Multi-shell Case Using Legendre Expansion
    %
    Rq_mag = rownorm(Rq);
    %
    if ~all( Rq_mag < R(1)+eps )    % check if dipoles within the brain
        warndlg('Legendre method assumes all dipoles(s) inside brain layer - please modify dipole locations OR use the "Berg" approach')
        return
    end
    %
    %  compute weights fn. fn depends only on the radii and cdv
    Rq_mag = rownorm(Rq);
    Re_mag = R(NL);         % Radius of outermost layer (Sensor distance from origin
    %
    %
    if nargin < 5,     % check # of inputs to see if nmax was specified
        nmax = fix(10/(1-max(Rq_mag)/Re_mag));    % Default for # Legendre Series Terms
    end
    %
    Ren = Re/Re_mag;
    Rqn = Rq./[Rq_mag,Rq_mag,Rq_mag];
    for k = 1:NL-1 
        s(k) = sigma(k)/sigma(k+1);
    end
    a = Re_mag./R;
    ainv = R/Re_mag;
    sm1 = s-1;
    twonp1 = 2*[1:nmax]+1;
    twonp1 = twonp1(:);
    f = zeros(nmax,1);
    %
    for n = 1:nmax
        np1 = n+1;
        Mc = eye(2);
        for k = 2:NL-1,
            Mc = Mc*[n+np1*s(k),  np1*sm1(k)*a(k)^twonp1(n);...
                    n*sm1(k)*ainv(k)^twonp1(n) , np1+n*s(k)];
        end;
        Mc(2,:) = [n*sm1(1)*ainv(1)^twonp1(n) , np1+n*s(1)]*Mc;
        Mc = Mc/(twonp1(n))^(NL-1);
        f(n) = n/(n*Mc(2,2)+np1*Mc(2,1));
    end;
    %
    onevec = ones(M,1);
    wtemp = ((twonp1./[1:nmax]').*f)/(4*pi*sigma(NL)*R(NL)^2); 
    n = [1:nmax]';
    nm1 = n-1;
    
    for i = 1:P,  % loop over all dipoles
        rqn = [Rqn(i,1)*onevec,Rqn(i,2)*onevec,Rqn(i,3)*onevec];
        cosgamma = dotprod(rqn,Ren); 
        rqn = rqn(1,:);             
        [Pl,dP] = dlegpoly(nmax,cosgamma);  % evaluate legendre poly and its derivative
        ratio = (Rq_mag(i)/Re_mag).^nm1;
        z = Ren- cosgamma*rqn;
        w = wtemp.*ratio;
        
        Gterm1 = Pl'*(w.*n);
        Gterm2 = dP'*w;
        G(:,3*i-2:3*i) = Gterm1*rqn + [z(:,1).*Gterm2,z(:,2).*Gterm2,z(:,3).*Gterm2];
    end
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%----------------------------------------------------------------------------------%%
    %%% This part computes the potential for a dipole contained within a multi-layer
    %%% isontropic sphere using Berg Parameter Approximation (Zhang Eqs 1i',5i" and 6i)
    %%% The EEG single-shell solution uses the vector form (which avoids computationally
    %%% expensive intrinsic functions) described by Mosher et al ("EEG and MEG: Forward 
    %%% Solutions for inverse problems" IEEE BME Trans March 1999)  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
elseif (NL>1)&(method==2)
    %
    if (~exist('mu_berg_in')|~exist('lam_berg_in'))
        error('Berg Parameters have not been specified!!!')
    elseif size(mu_berg_in)~=size(lam_berg_in)
        error('Berg Parameters are of Unequal Lengths!!!')
    else
        J = length(mu_berg_in);
        mu_berg = mu_berg_in;
        lam_berg = lam_berg_in;
    end
    %
    Re_mag = repmat(R(NL),P,M);           %(PxM)
    Re_mag_sq = repmat(R(NL)*R(NL),P,M);  %(PxM)
    Rq_mag = rownorm(Rq);                 %(Px1)
    Rq_mag_sq = Rq_mag.*Rq_mag;           %(Px1)
    Re_dot_Rq = Rq*Re';                   %(PxM)
    %
    for k=1:J
        %
        Rq1 = mu_berg(k)*Rq;                            %(Px3)
        Rq1_mag = mu_berg(k)*Rq_mag;                    %(Px1)
        Rq1_mag_sq = (mu_berg(k)*mu_berg(k))*Rq_mag_sq; %(Px1)
        Re_dot_Rq1 = mu_berg(k)*Re_dot_Rq;              %(PxM)
        %
        const = 4.0*pi*sigma(NL);
        const1 = const/lam_berg(k);
        term = 1./(const1*Rq1_mag_sq);                  %(PxM)
        %
        %%% This part checks for the presence of Berg dipoles which are external to
        %%% the sphere. For those dipoles external to the sphere, the dipole parameters
        %%% are replaced with the electrical image (internal to sphere) of the dipole
        %
        nx = find(Rq1_mag > R(NL));
        %
        if nx>0
            Rq1_temp = Rq1(nx,:);
            Rq1(nx,:) = R(NL)*R(NL)*Rq1_temp./repmat((rownorm(Rq1_temp).*rownorm(Rq1_temp)),1,3);
            Rq1_mag(nx,1) = rownorm(Rq1(nx,:));
            Rq1_mag_sq(nx,1) = Rq1_mag(nx,1).*Rq1_mag(nx,1);
            Re_dot_Rq1(nx,:) = R(NL)*R(NL)*Re_dot_Rq1(nx,:)./repmat((rownorm(Rq1_temp).*rownorm(Rq1_temp)),1,M);
            term(nx,:) = (R(NL)/rownorm(Rq1_temp))'.*term(nx,:);
        end
        %
        %%% Calculation of Forward Gain Matrix Contribution due to K-th Berg Dipole
        %
        Rq1_mag = repmat(Rq1_mag,1,M);                %(PxM)
        Rq1_mag_sq = repmat(Rq1_mag_sq,1,M);          %(PxM)
        term = repmat(term,1,M);                      %(PxM)
        %
        d_mag = reshape( rownorm(reshape(repmat(Re,1,P)',3,P*M)' ...
            -repmat(Rq1,M,1)) ,P,M);        %(PxM)
        d_mag_cub = d_mag.*d_mag.*d_mag;                        %(PxM)
        %
        F_scalar = d_mag.*(Re_mag.*d_mag+Re_mag_sq-Re_dot_Rq1); %(PxM)
        %
        c1 = term.*(2*( (Re_dot_Rq1-Rq1_mag_sq)./d_mag_cub) ...
            + 1./d_mag - 1./Re_mag);                 %(PxM)
        c2 = term.*((2./d_mag_cub) + (d_mag+Re_mag)./(Re_mag.*F_scalar));
        %
        G = G +    reshape(repmat((c1 - c2.*Re_dot_Rq1)',3,1),M,3*P) ...
            .*repmat(reshape(Rq1',1,3*P),M,1) ...
            +  reshape(repmat((c2.*Rq1_mag_sq)',3,1),M,3*P) ...
            .*repmat(Re,1,P);
    end
    %
    %---------------------------------------------------------------------------------------
    %
end    % End of Check for Forward Model Calculation Method
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
