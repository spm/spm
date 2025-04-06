function [D] = spm_eeg_invert_classic(D,val)
% Parallelized version of spm_eeg_invert_classic
% This version processes multiple time windows (wois) in parallel

Nl = length(D);

if Nl>1
    error('function only defined for a single subject');
end

if nargin > 1
    D.val = val;
elseif ~isfield(D, 'val')
    D.val = 1;
end

val=D.val;
inverse = D.inv{val}.inverse;

try, type = inverse.type;   catch, type = 'GS';     end
try, s    = inverse.smooth; catch, s    = 0.6;      end
try, Np   = inverse.Np;     catch, Np   = 256;      end
try, Nr   = inverse.Nr;     catch, Nr   = 16;       end
try, xyz  = inverse.xyz;    catch, xyz  = [0 0 0];  end
try, rad  = inverse.rad;    catch, rad  = 128;      end
try, hpf  = inverse.hpf;    catch, hpf  = 48;       end
try, lpf  = inverse.lpf;    catch, lpf  = 0;        end
try, sdv  = inverse.sdv;    catch, sdv  = 4;        end
try, Han  = inverse.Han;    catch, Han  = 1;        end
try, woi  = inverse.woi;    catch, woi  = [];       end
try, Nm   = inverse.Nm;     catch, Nm   = [];       end
try, Nt   = inverse.Nt;     catch, Nt   = [];       end
try, Ip   = inverse.Ip;     catch, Ip   = [];       end
try, QE    = inverse.QE;     catch,  QE=1;          end
try, Qe0   = inverse.Qe0;     catch, Qe0   = exp(-5);       end
try, inverse.A;     catch, inverse.A   = [];       end
try, SHUFFLELEADS=inverse.SHUFFLELEADS;catch, SHUFFLELEADS=0;end;

type = inverse.type;

modalities = D.inv{val}.forward.modality;

Nmax  = 16;

fprintf('Checking leadfields')
[L,D] = spm_eeg_lgainmat(D);
Nd=size(L,2);

if ~isempty(Ip)
    Np   = length(Ip);
else
    Ip=ceil([1:Np]*Nd/Np);
end

persistent permind;

rand(2)
if SHUFFLELEADS
    rng('shuffle')
    if isempty(permind)
        permind=randperm(size(L,1));
    end
    L=L(permind,:);
    warning('PERMUTING LEAD FIELDS !');
    permind(1:3)
end

if size(modalities,1)>1
    error('not defined for multiple modalities');
end
Ic  = setdiff(D.indchantype(modalities), badchannels(D));
Nd    = size(L,2);

fprintf(' - done\n')

if s>=1
    smoothtype='mesh_smooth';
else
    smoothtype='msp_smooth';
end
if s<0
    smoothtype='mm_smooth';
    s=-s;
end
vert  = D.inv{val}.mesh.tess_mni.vert;
face  = D.inv{val}.mesh.tess_mni.face;
M1.faces=face;
M1.vertices=vert;

switch smoothtype
    case 'mesh_smooth'
        fprintf('Using SPM smoothing for priors:')
        
        Qi    = speye(Nd,Nd);
        [QG]=spm_mesh_smooth(M1,Qi,round(s));
        QG    = QG.*(QG > exp(-8));
        
        QG    = QG*QG;
        disp('Normalising smoother');
        QG=QG./repmat(sum(QG,2),1,size(QG,1));
        
    case 'msp_smooth'
        fprintf('Computing Green function from graph Laplacian to smooth priors:')
        
        A     = spm_mesh_distmtx(struct('vertices',vert,'faces',face),0);
        
        GL    = A - spdiags(sum(A,2),0,Nd,Nd);
        GL    = GL*s/2;
        Qi    = speye(Nd,Nd);
        QG    = sparse(Nd,Nd);
        
        for i = 1:8
            QG = QG + Qi;
            Qi = Qi*GL/i;
        end
        
        QG    = QG.*(QG > exp(-8));
        QG    = QG*QG;
        
    case 'mm_smooth'
        
         kernelname=spm_eeg_smoothmesh_mm(D.inv{val}.mesh.tess_ctx,s);
         asmth=load(kernelname,'QG','M','faces');
         if isfield(asmth, 'M')
             if isa(asmth.M, 'gifti')
                 meshFaces = double(asmth.M.faces);
             else
                 error('M is not a recognized format.');
             end
         elseif isfield(asmth, 'faces')
             meshFaces = double(asmth.faces);
         else
             error('No valid mesh face data found in the smoothing kernel file.');
         end
         if max(int32(meshFaces)-int32(face))~=0
             error('Smoothing kernel used different mesh');
         end
         QG=asmth.QG;
         clear asmth;
          
end

clear Qi A GL
fprintf(' - done\n')

fprintf('Optimising and aligning spatial modes ...\n')

if isempty(inverse.A)
    if isempty(Nm)
        [U,ss,vv]    = spm_svd((L*L'),exp(-16));
        A     = U';
        UL    = A*L;
        
    else
        [U,ss,vv]    = spm_svd((L*L'),0);
        if length(ss)<Nm
            disp('number available');
            length(ss)
            error('Not this many spatial modes in lead fields');
        end
        
        ss=ss(1:Nm);
        disp('using preselected number spatial modes !');
        A     = U(:,1:Nm)';
        UL    = A*L;
    end
else
    disp('Using pre-specified spatial modes');
    if isempty(Nm)
        error('Need to specify number of spatial modes if U is prespecified');
    end
    A=inverse.A;
    UL=A*L;
end

Nm    = size(UL,1);

clear ss vv

fprintf('Using %d spatial modes',Nm)

Is    = 1:Nd;
Ns    = length(Is);

F=zeros(1,size(woi,1));

if isempty(woi)
    woi = 1000*[min(D.time) max(D.time)];
end
R2=zeros(1,size(woi,1));
VE=zeros(1,size(woi,1));

% Create a parallel pool if one doesn't exist
if isempty(gcp('nocreate'))
    parpool;
end

% Process each woi in parallel
parfor w_idx=1:size(woi,1)
    [F(w_idx), R2(w_idx), VE(w_idx), J_temp{w_idx}, M_temp{w_idx}, Cq_temp{w_idx}, ...
     U_temp{w_idx}, V_temp{w_idx}, Vq_temp{w_idx}, S_temp{w_idx}, It_temp{w_idx}, ...
     Ik_temp{w_idx}, ID_temp{w_idx}, pst_temp{w_idx}, dct_temp{w_idx}] = ...
        process_woi(D, val, w_idx, woi, A, UL, QG, Is, Ns, Ip, Np, QE, Qe0, type, Nm, Nmax, Nt, Nr, Han, lpf, hpf, sdv, Ic);
end

% Combine results from parallel processing
for w_idx=1:size(woi,1)
    if w_idx == 1
        inverse.type   = type;
        inverse.smooth = s;
        inverse.M      = M_temp{w_idx};
        inverse.J      = J_temp{w_idx};
        inverse.L      = UL;
        inverse.qC     = Cq_temp{w_idx};
        inverse.tempU  = U_temp{w_idx};
        inverse.V      = V_temp{w_idx};
        inverse.qV     = Vq_temp{w_idx};
        inverse.T      = S_temp{w_idx};
        inverse.U      = {A};
        inverse.Is     = Is;
        inverse.It     = It_temp{w_idx};
        inverse.Ik     = Ik_temp{w_idx};
        try
            inverse.Ic{1} = Ic;
        catch
            inverse.Ic = Ic;
        end
        inverse.Nd     = Nd;
        inverse.pst    = pst_temp{w_idx};
        inverse.dct    = dct_temp{w_idx};
        inverse.ID     = ID_temp{w_idx};
    end
    inverse.F(w_idx)   = F(w_idx);
    inverse.R2(w_idx)  = R2(w_idx);
    inverse.VE(w_idx)  = R2(w_idx).*VE(w_idx);
end

inverse.woi    = woi;
inverse.Ip     = Ip;
inverse.modality = modalities;

D.inv{val}.inverse = inverse;
D.inv{val}.method  = 'Imaging';

if ~spm('CmdLine')
    spm_eeg_invert_display(D);
    drawnow
end

end

function [F_out, R2_out, VE_out, J_out, M_out, Cq_out, U_out, V_out, Vq_out, ...
          S_out, It_out, Ik_out, ID_out, pst_out, dct_out] = ...
    process_woi(D, val, w_idx, woi, A, UL, QG, Is, Ns, Ip, Np, QE, Qe0, type, Nm, Nmax, Nt, Nr, Han, lpf, hpf, sdv, Ic)

% This function processes a single time window of interest (woi)
w = woi(w_idx,:);

if ~isempty(Ip)
    Np = length(Ip);
else
    Ip = ceil([1:Np]*Ns/Np);
end

AY    = {};
AYYA  = 0;

It = (w/1000 - D.timeonset)*D.fsample + 1;
It = max(1,It(1)):min(It(end), length(D.time));
It = fix(It);
disp(sprintf('Number of samples %d',length(It)));

pst = 1000*D.time;
pst = pst(It);
dur = (pst(end) - pst(1))/1000;
dct = (It - It(1))/2/dur;
Nb  = length(It);

K   = exp(-(pst - pst(1)).^2/(2*sdv^2));
K   = toeplitz(K);
qV  = sparse(K*K');

T   = spm_dctmtx(Nb,Nb);

j   = find((dct >= lpf) & (dct <= hpf));
T   = T(:,j);
dct = dct(j);

if Han
    W = sparse(1:Nb,1:Nb,spm_hanning(Nb));
else
    W = 1;
end

try
    trial = D.inv{D.val}.inverse.trials;
catch
    trial = D.condlist;
end
Ntrialtypes = length(trial);

YY = 0;
N  = 0;

badtrialind = D.badtrials;
Ik = [];
for j = 1:Ntrialtypes
    c = D.indtrial(trial{j});
    [c1,ib] = intersect(c,badtrialind);
    c = c(setxor(1:length(c),ib));
    Ik = [Ik c];
    Nk = length(c);
    for k = 1:Nk
        Y = A*D(Ic,It,c(k));
        YY = YY + Y'*Y;
        N = N + 1;
    end
end
YY = YY./N;

YY = W'*YY*W;
YTY = T'*YY*T;

if isempty(Nt)
    [U, E] = spm_svd(YTY,exp(-8));
    if isempty(U)
        warning('nothing found using spm svd, using svd');
        [U E] = svd(YTY);
    end
    E = diag(E)/trace(YTY);
    Nr = min(length(E),Nmax);
    Nr = max(Nr,1);
else
    [U, E] = svd(YTY);
    E = diag(E)/trace(YTY);
    disp('Fixed number of temporal modes');
    Nr = Nt;
end

V = U(:,1:Nr);
VE_out = sum(E(1:Nr));

fprintf('Using %i temporal modes, ',Nr)
fprintf('accounting for %0.2f percent average variance\n',full(100*VE_out))

S = T*V;
Vq = S*pinv(S'*qV*S)*S';

UYYU = 0;
AYYA = 0;
Nn = 0;
AY = {};
Ntrials = 0;

for j = 1:Ntrialtypes
    UY{j} = sparse(0);
    c = D.indtrial(trial{j});
    [c1,ib] = intersect(c,badtrialind);
    c = c(setxor(1:length(c),ib));
    Nk = length(c);
    
    for k = 1:Nk
        Y = D(Ic,It,c(k))*S;
        Y = A*Y;
        
        Nn = Nn + Nr;
        
        YY = Y*Y';
        Ntrials = Ntrials+1;
        
        UY{j} = UY{j} + Y;
        UYYU = UYYU + YY;
        
        AY{end + 1} = Y;
        AYYA = AYYA + YY;
    end
end

AY = spm_cat(AY);

ID = spm_data_id(AY);

AQeA = A*QE*A';
Qe{1} = AQeA/(trace(AQeA));

Q0 = Qe0*trace(AYYA)*Qe{1}./sum(Nn);

allind = [];
switch(type)
    case {'MSP','GS','ARD'}
        Qp = {};
        LQpL = {};
        for i = 1:Np
            q = QG(:,Ip(i));
            Qp{end + 1}.q = q;
            LQpL{end + 1}.q = UL*q;
        end
        
    case {'EBB'}
        disp('NB smooth EBB algorithm !');
        InvCov = spm_inv(AYYA);
        allsource = sparse(Ns,1);
        Sourcepower = sparse(Ns,1);
        for bk = 1:Ns
            q = QG(:,bk);
            
            smthlead = UL*q;
            if ~all(smthlead==0)
                normpower = 1/(smthlead'*smthlead);
                Sourcepower(bk) = 1/(smthlead'*InvCov*smthlead);
                allsource(bk) = Sourcepower(bk)./normpower;
            end
        end
        allsource = allsource/max(allsource);
        
        Qp{1} = diag(allsource);
        LQpL{1} = UL*diag(allsource)*UL';
        
    case {'EBBcorr'}
        disp('NB smooth correlated source EBB algorithm !');
        InvCov = spm_inv(AYYA);
        halfNs = ceil(Ns/2);
        allsource = sparse(Ns,1);
        Sourcepower = sparse(Ns,1);
        DualSourcepower = sparse(halfNs,1);
        alldualsource = sparse(Ns,1);
        for bk = 1:Ns
            q = QG(:,bk);
            
            smthlead = UL*q;
            normpower = 1/(smthlead'*smthlead);
            Sourcepower(bk) = 1/(smthlead'*InvCov*smthlead);
            allsource(bk) = Sourcepower(bk)./normpower;
        end
        
        leftbrainind = find(vert(:,1)<0);
        
        for bk = 1:length(vert)
            vertind = bk;
            
            reflectpos = [-vert(vertind,1) vert(vertind,2) vert(vertind,3)];
            d1 = vert - repmat(reflectpos,length(vert),1);
            dist1 = dot(d1',d1');
            [d1,reflectind] = min(dist1);
            q = QG(:,vertind)+QG(:,reflectind);
            smthlead = UL*q;
            normpower = 1/(smthlead'*smthlead);
            DualSourcepower(vertind) = 1/(smthlead'*InvCov*smthlead);
            alldualsource(vertind) = alldualsource(vertind)+DualSourcepower(vertind)./(normpower*4);
            alldualsource(reflectind) = alldualsource(reflectind)+DualSourcepower(vertind)./(normpower*4);
        end
        
        allsource = allsource+alldualsource;
        allsource = allsource/max(allsource);
        
        Qp{1} = diag(allsource);
        LQpL{1} = UL*diag(allsource)*UL';
        
    case {'EBBgs'}
        allsource = zeros(Ntrials,Ns);
        for ii = 1:Ntrials
            InvCov = spm_inv(YYep{ii});
            Sourcepower = zeros(Ns,1);
            for bk = 1:Ns
                normpower = 1/(UL(:,bk)'*UL(:,bk));
                Sourcepower(bk) = 1/(UL(:,bk)'*InvCov*UL(:,bk));
                allsource(ii,bk) = Sourcepower(bk)./normpower;
            end
            
            Qp{ii}.q = allsource(ii,:);
        end
        
    case {'LOR','COH'}
        Qp{1} = speye(Ns,Ns);
        LQpL{1} = UL*UL';
        
        Qp{2} = QG;
        LQpL{2} = UL*Qp{2}*UL';
        
    case {'IID','MMN'}
        Qp{1} = speye(Ns,Ns);
        LQpL{1} = UL*UL';
end

fprintf('Using %d spatial source priors provided\n',length(Qp));

QP = {};
LQP = {};
LQPL = {};

switch(type)
    case {'MSP','GS','EBBgs'}
        Np = length(Qp);
        Q = zeros(Ns,Np);
        for i = 1:Np
            Q(:,i) = Qp{i}.q;
        end
        Q = sparse(Q);
        
        MVB = spm_mvb(AY,UL,[],Q,Qe,16);
        
        Qcp = Q*MVB.cp;
        QP{end + 1} = sum(Qcp.*Q,2);
        LQP{end + 1} = (UL*Qcp)*Q';
        LQPL{end + 1} = LQP{end}*UL';
end

switch(type)
    case {'MSP','ARD'}
        [Cy,h,Ph,F_out] = spm_sp_reml(AYYA,[],[Qe LQpL],Nn);
        
        Ne = length(Qe);
        Np = length(Qp);
        
        hp = h(Ne + (1:Np));
        
        qp = sparse(0);
        for i = 1:Np
            if hp(i) > max(hp)/128
                qp = qp + hp(i)*Qp{i}.q*Qp{i}.q';
            end
        end
        
        QP{end + 1} = diag(qp);
        LQP{end + 1} = UL*qp;
        LQPL{end + 1} = LQP{end}*UL';
end

switch(type)
    case {'IID','MMN','LOR','COH','EBB','EBBcorr'}
        [Cy,h,Ph,F_out] = spm_reml_sc(AYYA,[],[Qe LQpL],Nn,-4,16,Q0);
        
        Ne = length(Qe);
        Np = length(Qp);
        
        hp = h(Ne + (1:Np));
        qp = sparse(0);
        for i = 1:Np
            qp = qp + hp(i)*Qp{i};
        end
        
        QP{end + 1} = diag(qp);
        LQP{end + 1} = UL*qp;
        LQPL{end + 1} = LQP{end}*UL';
end

fprintf('Inverting subject 1\n')

Np = length(LQPL);
Ne = length(Qe);

Q = [{Q0} LQPL];

if rank(AYYA)~=size(A,1)
    rank(AYYA);
    size(AYYA,1);
    warning('AYYA IS RANK DEFICIENT');
end

[Cy,h,Ph,F_out] = spm_reml_sc(AYYA,[],Q,Nn,-4,16,Q0);

Cp = sparse(0);
LCp = sparse(0);
hp = h(Ne + (1:Np));
for j = 1:Np
    Cp = Cp + hp(j)*QP{j};
    LCp = LCp + hp(j)*LQP{j};
end

M = LCp'/Cy;

Cq = Cp - sum(LCp.*M')';

SSR = 0;
SST = 0;
J = {};

for j = 1:Ntrialtypes
    J{j} = M*UY{j};
    
    SSR = SSR + sum(var((UY{j} - UL*J{j})));
    SST = SST + sum(var(UY{j}));
end

R2_out = 100*(SST - SSR)/SST;
fprintf('Percent variance explained %.2f (%.2f)\n',full(R2_out),full(R2_out*VE_out));

% Outputs for combining later
J_out = J;
M_out = M;
Cq_out = Cq;
U_out = U;
V_out = V;
Vq_out = Vq;
S_out = S;
It_out = It;
Ik_out = Ik;
ID_out = ID;
pst_out = pst;
dct_out = dct;
end