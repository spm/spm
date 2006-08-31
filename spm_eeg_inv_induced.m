function D = spm_eeg_inv_induced(D,Qe,Qp)

%=======================================================================
% Inverse Solution for induced (and evoked) power
%
% FORMAT D = spm_eeg_inv_induced(D,Qe,Qp)
% Input:
% D                 - input data structure
% Qe        - noise covariance structure
% Qp        - prior source covaraince structure
% Output:
% D                     - same data struct including the inverse solution files and variables
%=======================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout
% $Id: spm_eeg_inv_induced.m 607 2006-08-31 12:29:39Z james $


if length(D.events.code) ~= D.Nevents
   error(sprintf('Single trial data required\n'));
end
val = length(D.inv);

qe = length(Qe);
qp = length(Qp);

% LOAD GAIN MATRIX
variabl = load(D.inv{val}.forward.gainmat);
name    = fieldnames(variabl);
G       = getfield(variabl , name{1});


% SOURCE SPACE DIMENSION REDUCTION
if (D.inv{val}.mesh.Ctx_Nv - D.inv{val}.inverse.dim)
   Msize     = D.inv{val}.inverse.dim
   if isempty(D.inv{val}.inverse.priors.level2{3}.filename)
       error(sprintf('Multivariate Source Prelocalisation has not been performed\n'));
   end
   variabl   = load(D.inv{val}.inverse.priors.level2{3}.filename);
   name      = fieldnames(variabl);
   APM       = getfield(variabl , name(1));
   [APMs,Is] = sort(-APM);
   APMs      = -APMs;
   Is        = Is(1:Msize);
   G         = G(:,Is);
   GQpG      = {};
   for i = 1:length(Qp)
       Qp{i}   = Qp{i}(Is,Is);
   end
end


% VARIANCE COMPONENTS
for i = 1:length(Qp)
   GQpG{i} = G*Qp{i}*G';
   GQpG{i} = GQpG{i}/norm(GQpG{i},'fro');
end

for i = 1:length(Qe)
   Qe{i}   = Qe{i}/norm(Qe{i},'fro');
end

woi        = D.inv{val}.inverse.woi;
contrast   = D.inv{val}.inverse.contrast;
Nsens      = D.Nchannels;
Nsour      = size(G,2);

woi(1) = round(woi(1)*(D.Radc/1000)) + D.events.start + 1;
woi(2) = round(woi(2)*(D.Radc/1000)) + D.events.start + 1;
Nsamp  = woi(2) - woi(1) + 1;

Ic = find(contrast);
if length(Ic) > 1
   error('You should select one trial type onlym\n');
end
It = find(D.events.code == D.events.types(Ic));
Ntrial = length(It);
Yi = [];
for i = 1:Ntrial
Yi = [Yi squeeze(D.data(:,woi(1):woi(2),It(i)))];
end


% SOURCE TEMPORAL STRUCTURE
V     = eye(Nsamp,Nsamp);        % temporal correlations
T     = convmtx(spm_Npdf(-8:8,0,2^2),Nsamp);
T     = T(:,[1:Nsamp] + 8);                           % temporal dispersion
w     = diag([1:Nsamp].^(9/8).*exp(-[1:Nsamp]*3.1641/Nsamp));     % window
[S v] = eig(w*T*T'*w);
dv    = diag(v)/sum(diag(v));
dr = 0; ir = 0;
while dr < 0.95
   ir = ir + 1;
   dr = dr + dv(ir);
end
r     = ir;                      % dimension of signal subspace
clear dv dr ir
S     = S(:,[1:r]);
SVS   = S'*V*S;                                       % correlation (signal)

   % time-frequency subspace
Fband = D.inv{val}.inverse.fboi(2) - D.inv{val}.inverse.fboi(1) + 1;
Fstep = floor((Fband-1)/10) + 1;
WT = 2*pi*[1:Nsamp]'/D.Radc;
W = [];
Fc = D.inv{val}.inverse.fboi(1);
while Fc < D.inv{val}.inverse.fboi(2)
   W = [W sin(Fc*WT) cos(Fc*WT)];
   Fc = Fc + Fstep;
end
v = spm_Npdf(1:Nsamp,round(Nsamp/2),(Nsamp/4)^2); % time window
W = diag(v)*W;
clear WT Fc Fstep Fband


% EVOKED RESPONSE
Ye = Yi*kron(ones(Ntrial,1),speye(Nsamp))/Ntrial;
if D.inv{val}.inverse.activity == 'evoked & induced';
   YY = Ye*S*inv(SVS)*S'*Ye'/r;
   X  = ones(Nsens,1);
   Q  = {Qe{:} GQpG{:}};

   % ReML hyperparameter estimates
   [Cev,hev,Phev,Fev] = spm_reml(YY,X,Q,Nsamp,1);
   clear YY Q

   % Display ReML outcome
   disp(sprintf('\n'));
   disp('--- EVOKED RESPONSE ---');
   disp(['Model Log-Evidence:  ' num2str(Fev)]);
   DisplayString = ['Noise covariance hyperparameters:  '];
   for i = 1:qe
       DisplayString = [DisplayString '  ' num2str(hev(i))];
   end
   disp(DisplayString);
   DisplayString = ['Source covariance hyperparameters:  '];
   for i = 1:qp
       DisplayString = [DisplayString '  ' num2str(hev(qe + i))];
   end
   disp(DisplayString);

   % MAP parameter estimates
   Ce = sparse(Nsens,Nsens);
   for i = 1:qe
       Ce = Ce + hev(i)*Qe{i};
   end
   Cp = sparse(Nsour,Nsour);
   for i = 1:qp
       Cp = Cp + hev(qe + i)*Qp{i};
   end

   CpG  = Cp*G';
   GCpG = G*CpG;
   MAP  = CpG*inv(GCpG + Ce);
   C    = inv(G'*inv(Ce)*G + inv(Cp + speye(Nsour,Nsour)*1e-6));
   clear CpG GCpG Ce Cp

       % instantaneous source activity
   J = MAP*Ye*S*S';
   if (D.inv{val}.mesh.Ctx_Nv - D.inv{val}.inverse.dim)
       Jfull       = sparse(Nsour,Nsamp);
       Jfull(Is,:) = J;
       Jev         = Jfull;
   else
       Jev         = J;
   end
   clear J Jfull

       % cross-energy in channel space
   K    = S*S'*W*W'*S*S';
   Eev  = Ye*K*Ye';
   Eev  = diag(Eev);

       % cross-energy in source space
   Gev  = MAP*Eev*MAP' + C*trace(K*V);
   Gev  = diag(Gev);
end


% INDUCED RESPONSE
Yi = Yi - kron(ones(1,Ntrial),Ye);
clear Ye
YY = Yi*kron(speye(Ntrial),S*inv(SVS)*S')*Yi'/(Ntrial*r);
X  = ones(Nsens,1);
Q  = {Qe{:} GQpG{:}};

% ReML hyperparameter estimates
[Cind,hind,Phind,Find] = spm_reml(YY,X,Q,Nsamp,1);
clear YY Q

% Display ReML outcome
disp(sprintf('\n'));
disp('--- INDUCED RESPONSE ---');
disp(['Model Log-Evidence:  ' num2str(Find)]);
DisplayString = ['Noise covariance hyperparameters:  '];
for i = 1:qe
   DisplayString = [DisplayString '  ' num2str(hind(i))];
end
disp(DisplayString);
DisplayString = ['Source covariance hyperparameters:  '];
for i = 1:qp
   DisplayString = [DisplayString '  ' num2str(hind(qe + i))];
end
disp(DisplayString);

% MAP estimate of the energy
Ce = sparse(Nsens,Nsens);
for i = 1:qe
   Ce = Ce + hind(i)*Qe{i};
end
Cp = sparse(Nsour,Nsour);
for i = 1:qp
   Cp = Cp + hind(qe + i)*Qp{i};
end
clear Qe Qp

CpG     = Cp*G';
GCpG    = G*CpG;
MAP     = CpG*inv(GCpG + Ce);
ToBeInv = G'*inv(Ce)*G;
clear G
Cp      = Cp + 1e-6*speye(Nsour,Nsour);
Cp      = inv(Cp);
ToBeInv = ToBeInv + Cp;
clear Cp
C       = inv(ToBeInv);
% C     = inv(G'*inv(Ce)*G + inv(Cp + speye(Nsour,Nsour)*1e-6));
SSVSS = S*SVS*S';
clear CpG GCpG Ce ToBeInv


   % cross-energy in channel space
K    = S*S'*W*W'*S*S';
Eind = Yi*kron(speye(Ntrial),K)*Yi'/Ntrial;
Eind = diag(Eind);
clear Yi Ce Cp YY X

   % cross-energy in source space
Gind = MAP*Eind*MAP' + C*trace(K*V);
Gind = diag(Gind);


% Save results
if strcmp(D.inv{val}.inverse.activity,'induced')

   [pth,nam,ext] = spm_fileparts(D.fname);
   woi           = D.inv{val}.inverse.woi;
   Ntime         = clock;
   D.inv{val}.inverse.resfile = [nam '_remlmat_' num2str(woi(1)) '_' num2str(woi(2)) 'ms_induced_' num2str(Ntime(4)) 'H' num2str(Ntime(5)) '.mat'];
   D.inv{val}.inverse.LogEv   = Find;
   if str2num(version('-release'))>=14
       save(fullfile(pth,D.inv{val}.inverse.resfile), '-V6','Cind','hind','Phind','Find','Eind','Gind');
   else
       save(fullfile(pth,D.inv{val}.inverse.resfile),'Cind','hind','Phind','Find','Eind','Gind');
   end
   clear Cind hind Phind Find Eind Gind

   % Temporary Visualization
   Rind = diag(Gind);
   load(D.inv{val}.mesh.tess_ctx);
   colormap jet
   axis off
   patch('Vertices',vert,'Faces',face,'FaceVertexCData',Rind,'FaceColor','flat');
   shading interp
   colorbar
   title('Induced Power');
   view(-90,0);
   cameramenu

else

   [pth,nam,ext] = spm_fileparts(D.fname);
   woi           = D.inv{val}.inverse.woi;
   Ntime         = clock;
   D.inv{val}.inverse.LogEv = [Fev Find];
   D.inv{val}.inverse.resfile{1} = [nam '_remlmat_' num2str(woi(1)) '_' num2str(woi(2)) 'ms_evoked' num2str(Ntime(4)) 'H' num2str(Ntime(5)) '.mat'];
   if str2num(version('-release'))>=14
       save(fullfile(pth,D.inv{val}.inverse.resfile{1}), '-V6','Cev','hev','Phev','Fev','Jev','Eev','Gev');
   else
       save(fullfile(pth,D.inv{val}.inverse.resfile{1}),'Cev','hev','Phev','Fev','Jev','Eev','Gev');
   end
   clear Cev hev Phev Fev Jev Eev Gev

   % Temporary Visualization
   Rev = diag(Gev);
   load(D.inv{val}.mesh.tess_ctx);
   colormap jet
   axis off
   patch('Vertices',vert,'Faces',face,'FaceVertexCData',Rev,'FaceColor','flat');
   view(-90,0);
   shading interp
   colorbar
   title('Evoked Power');
   cameramenu


   D.inv{val}.inverse.resfile{2} = [nam '_remlmat_' num2str(woi(1)) '_' num2str(woi(2)) 'ms_induced' num2str(Ntime(4)) 'H' num2str(Ntime(5)) '.mat'];
   if str2num(version('-release'))>=14
       save(fullfile(pth,D.inv{val}.inverse.resfile{2}), '-V6','Cind','hind','Phind','Find','Eind','Gind');
   else
       save(fullfile(pth,D.inv{val}.inverse.resfile{2}),'Cind','hind','Phind','Find','Eind','Gind');
   end
   clear Cind hind Phind Find Eind Gind

   % Temporary Visualization
   spm_figure
   Rind = diag(Gind);
   axis off
   patch('Vertices',vert,'Faces',face,'FaceVertexCData',Rind,'FaceColor','flat');
   view(-90,0);
   shading interp
   colorbar
   title('Induced Power');
   cameramenu

end

if str2num(version('-release'))>=14
   save(fullfile(D.path,D.fname), '-V6','D');
else
   save(fullfile(D.path,D.fname),'D');
end