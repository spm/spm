function P = spm_eeg_inv_vbecd(P)
% Model inversion routine for ECDs using "variational Bayes"
%
% FORMAT P = spm_eeg_inv_vbecd(P)
%
% Input:
% structure P with fields:
%  forward      - structure containing the forward model, i.e. the "vol"
%                 and "sens" structure in a FT compatible format
%  bad          - list of bad channels, not to use.
%  y            - data vector
%  Nc           -
%  Niter        - maximum number of iterations
%  threshold_dF - threshold on free energy improvement to stop iterating
%  priors       - priors on parameters, hard and soft, as filled in (and
%                 described) in spm_eeg_inv_vbecd_gui.m.
%
% Output:
% same structure with extra fields
%  init         - initial valuse used for mu_w/s
%  ok           - flags indicating if everything was ok
%  dF           - successive (relative) improvement of F
%  post         - posterior value of estimated parameters and ther variance
%  Fi           - successive values of F
%  F            - Free energy final value.
%
% Reference:
% Kiebel et al., Variational Bayesian inversion of the equivalent current
% dipole model in EEG/MEG., NeuroImage, 39:728-741, 2008
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Christophe Phillips & Stefan Kiebel
% $Id: spm_eeg_inv_vbecd.m 3034 2009-04-01 15:12:55Z jean $




% unpack model, priors, data
%---------------------------

a10 = P.priors.a10; b10 = P.priors.b10;
a20 = P.priors.a20; b20 = P.priors.b20;
a30 = P.priors.a30; b30 = P.priors.b30;

Nd = length(P.priors.mu_w0)/3;
mu_w0 = P.priors.mu_w0;
mu_s0 = P.priors.mu_s0;
iS_w0 = P.priors.iS_w0;
iS_s0 = P.priors.iS_s0;

try
    Tw = P.priors.Tw;
catch
    Tw = eye(length(mu_w0));
end

try
    Ts = P.priors.Ts;
catch
    Ts = eye(length(mu_s0));
end

Vw = full(spm_svd(Tw));
Vs = full(spm_svd(Ts));

y = P.y;
y = y-mean(y); % Re-reference data to mean

P.dv = 10^-2; % used to compute step-size for gradients


%---------------
% initialization
%---------------

% ensure data have apropriate scale
% if strcmp(P.modality,'MEG')
%     % sc_y = norm(y)/10;
%     sc_y = 1e-8;
%     y = y/sc_y;
% else
%     sc_y = 1;
% end
sc_y = 1;
% Initialize posterior with prior
b10 = b10*sc_y.^2;
a1 = a10;
b1 = b10;
a2 = a20;
b2 = b20;
a3 = a30;
b3 = b30;
[u,s,v] = svd((a2/b2)*iS_w0);
mu_w = mu_w0 + 1e-8*u*diag(sqrt(diag(s)))*v'*randn(size(mu_w0));
S_w = u*diag(diag((s+eps).^-1))*v';
[u,s,v] = svd((a3/b3)*iS_s0);
mu_s = mu_s0 + 1e-8*u*diag(sqrt(diag(s)))*v'*randn(size(mu_s0));
S_s = u*diag(diag((s+eps).^-1))*v';

% get lead fields
[gmn, gm, dgm] = spm_eeg_inv_vbecd_getLF(mu_s, P.forward.sens, P.forward.vol,...
    P.dv.*ones(1, length(mu_w0)));

[Nc, Np] = size(gmn);
DE = kron(S_w+mu_w*mu_w', eye(Nc));
P.init.mu_s = mu_s;
P.init.mu_w = mu_w;


% Calulate free energy with priors
F(1) = -Nc/2*log(2*pi) + Nc/2*(psi(a1) - log(b1))...
    -a1/(2*b1)*(y'*y - 2*mu_w'*gmn'*y + gm'*DE*gm + trace(S_s*dgm'*DE*dgm))...
    -spm_kl_normal(Vw'*mu_w, Vw'*S_w*Vw, Vw'*mu_w0, Vw'*b2/a2*inv(iS_w0)*Vw)...
    -spm_kl_normal(Vs'*mu_s, Vs'*S_s*Vs, Vs'*mu_s0, Vs'*b3/a3*inv(iS_s0)*Vs)...
    -spm_kl_gamma(1/b1,a1,1/b10,a10)...
    -spm_kl_gamma(1/b2,a2,1/b20,a20)...
    -spm_kl_gamma(1/b3,a3,1/b30,a30);

P.gmn = gmn;


P.handles.hfig  = spm_figure('GetWin','Graphics');
spm_clf(P.handles.hfig)
P.handles.SPMdefaults.col = get(P.handles.hfig,'colormap');
P.handles.SPMdefaults.renderer = get(P.handles.hfig,'renderer');
set(P.handles.hfig,'userdata',P)

[P] = displayVBupdate(y,a1,b1,a2,b2,a3,b3,mu_w,mu_s,S_s,S_w,P,1,[],F);

P.ok = 1;
for i = 1:P.Niter
    
    % orientation parameters w
    SD = zeros(Np, Np);
    Dm = dgm*S_s*dgm';
    for j = 1:Nc
        ind = j + [0:Np-1]*Nc;
        d = Dm(ind, ind);
        SD = SD + d;
    end
    S_w = Tw*pinv(a2/b2*iS_w0 + a1/b1*(gmn'*gmn + SD))*Tw';
    mu_w = Tw*S_w*(a1/b1*gmn'*y + a2/b2*iS_w0*mu_w0);

    % location parameters s
    [mu_s,S_s,DE,gmn,gm,dgm,P] = ...
        modifiedGN4s(mu_s,S_s,iS_s0,mu_s0,mu_w,S_w,a1,b1,a3,b3,y,P,Ts);

    % precision on y
    a1 = Nc/2 + a10;
    b1 = 0.5*(y'*y...
        - 2*mu_w'*gmn'*y...
        + gm'*DE*gm...
        + trace(S_s*dgm'*DE*dgm))...
        + b10;

    % precision on w
    a2 = size(Vw,2)/2 + a20;
    b2 = 0.5*((Vw'*(mu_w-mu_w0))'*Vw'*iS_w0*Vw*(Vw'*(mu_w-mu_w0)) + ...
        trace(Vw'*iS_w0*S_w*Vw)) + b20;

    % precision on s
    a3 = size(Vs,2)/2 + a30;
    b3 = 0.5*((Vs'*(mu_s-mu_s0))'*Vs'*iS_s0*Vs*(Vs'*(mu_s-mu_s0)) + ...
        trace(Vs'*iS_s0*Vs*Vs'*S_s*Vs)) + b30;

    % compute neg free energy
    F(i+1) = -Nc/2*log(2*pi) + Nc/2*(psi(a1) - log(b1))...
        -a1/(2*b1)*(y'*y - 2*mu_w'*gmn'*y + gm'*DE*gm + trace(S_s*dgm'*DE*dgm))...
        -spm_kl_normal(Vw'*mu_w, Vw'*S_w*Vw, Vw'*mu_w0, Vw'*b2/a2*inv(iS_w0)*Vw)...
        -spm_kl_normal(Vs'*mu_s, Vs'*S_s*Vs, Vs'*mu_s0, Vs'*b3/a3*inv(iS_s0)*Vs)...
        -spm_kl_gamma(1/b1,a1,1/b10,a10)...
        -spm_kl_gamma(1/b2,a2,1/b20,a20)...
        -spm_kl_gamma(1/b3,a3,1/b30,a30);

    % update display
    P.gmn = gmn;
    [P] = displayVBupdate(y,a1,b1,a2,b2,a3,b3,mu_w,mu_s,S_s,S_w,P,i+1,[],F(end));
    pause(0.1)

    % check termination condition
    dF = (F(i+1)-F(i));%/abs(F(i));
    if abs(dF) < P.threshold_dF
        P.dF(i) = dF;
        str = sprintf('%3d/%d, F: %f\t dFr: %f', i, P.Niter, F(i+1), dF);
        fprintf('%s\n', str)
        break;
    end
    P.dF(i) = dF;
    str = sprintf('%3d/%d, F: %f\t dFr: %f', i, P.Niter, F(i+1), dF);
    fprintf('%s\n', str)

end

try
    set(P.handles.hte(2),'string','VB for ECDs: done.')
catch
    P.ok = 0;
end
set(P.handles.hfig,'userdata',[])
P = rmfield(P,'handles');


% rescale back to original units
mu_w = mu_w*sc_y;
S_w = S_w*sc_y^2;

% save results
post.mu_w = mu_w;
post.mu_s = mu_s;
post.S_w = S_w;
post.S_s = S_s;
post.a1 = a1; post.b1 = b1;
post.a2 = a2; post.b2 = b2;
post.a3 = a3; post.b3 = b3;
post.gmn = gmn; % lead field

P.post = post;
P.Fi = F;
P.F = F(end);


function [mu,Sigma,DE,gmn,gm,dgm,P] = ...
    modifiedGN4s(mu,S_s,iS_s0,mu_s0,mu_w,S_w,a1,b1,a3,b3,y,P,Ts)


% Compute variational energy and GN step from previous mode
PreviousMu = mu;
[PreviousI,Sigma,deltaMu,DE,gmn,gm,dgm] = ...
    logVarQ(PreviousMu,S_s,iS_s0,mu_s0,mu_w,S_w,a1,b1,a3,b3,y,P,Ts);

maxIter = 16;
rdI = 1e-4;
stop = 0;
it = 0;
while stop == 0
    it = it+1;
    % make a move
    mu = Ts*(PreviousMu + deltaMu);
    % get variational energy (as well as next move)
    [I,nextSigma,NextdeltaMu,nextDE,nextgmn,nextgm,nextdgm] = ...
        logVarQ(mu,S_s,iS_s0,mu_s0,mu_w,S_w,a1,b1,a3,b3,y,P,Ts);
    % calculate relative variational energy improvement
    deltaI = I-PreviousI;
    % check whether to stop, to accept move or to halve step
    if it <= maxIter && abs(deltaI./PreviousI)>rdI
        if deltaI<0     % halve step size
            deltaMu = 0.5*deltaMu;
            P.pc = 100*it./maxIter;
            P = displayVBupdate(y,a1,b1,[],[],a3,b3,mu_w,mu,Sigma,S_w,P,[],'mGN');
        else            % accept move
            % propose a new GN move from there on
            deltaMu = NextdeltaMu;
            DE = nextDE;
            gmn = nextgmn;
            gm = nextgm;
            dgm = nextdgm;
            Sigma = Ts*nextSigma*Ts';
            PreviousMu = mu;
            PreviousI = I;
            P.pc = 100*it./maxIter;
            P.gmn = gmn;
            P = displayVBupdate(y,a1,b1,[],[],a3,b3,mu_w,mu,Sigma,S_w,P,[],'ecd');
        end
    else                % stop Gauss-Newton search
        stop = 1;
    end
end
mu = PreviousMu;    % ensures mode is not updated if no improvement

function [I,Sigma,deltaMu,DE,gmn,gm,dgm] = ...
    logVarQ(PreviousMu,S_s,iS_s0,mu_s0,mu_w,S_w,a1,b1,a3,b3,y,P,Ts)
            
Nc = size(y,1);
In = speye(Nc);

vol = P.forward.vol;
sens = P.forward.sens;
dv = P.dv;

% first check if inside head
outside = ~forwinv_inside_vol(reshape(PreviousMu,3,[])',vol);
if ~all(~outside)
    I = -Inf;
    deltaMu = mu_s0-PreviousMu+eps;  % to correct for unluky initialization
    Sigma = [];
    DE = [];
    gmn = [];
    gm = [];
    dgm = [];
else
    % Compute lead fields (and gradients) at the mode
    [gmn, gm, dgm] =...
        spm_eeg_inv_vbecd_getLF(PreviousMu, sens, vol,dv.*sqrt(diag(S_s)));
    % Compute variational energy at the mode
    Gw = gmn*mu_w;
    dmu0 = mu_s0-PreviousMu;
    iSdmu0 = iS_s0*dmu0;
    dy = y-Gw;
    I = -0.5*(a3/b3)*dmu0'*iSdmu0 ...
        -0.5*(a1/b1)*( dy'*dy + trace(gmn'*gmn*S_w) );
    % Get local curvature of variational energy (-> cov matrix)
    DE = kron(S_w+mu_w*mu_w', eye(Nc));
    Sigma = pinv(a3/b3*iS_s0 + a1/b1*(dgm'*DE*dgm));
    % standard Gauss-Newton move from mode and curvature
    deltaMu = Sigma*( (a3/b3)*iSdmu0 ...
        + (a1/b1)*dgm'*(kron(mu_w,dy)) );%+ kron(S_w,In)*gm) );
end




function [P] = displayVBupdate(y,a1,b1,a2,b2,a3,b3,mu_w,mu_s,S_s,S_w,P,it,flag,F)

if ~exist('flag','var')
    flag = [];
end


if isempty(flag) || isequal(flag,'ecd')
    % plot dipoles
    try
        opt.ParentAxes = P.handles.axesECD;
        opt.hfig = P.handles.hfig;
        opt.handles.hp = P.handles.hp;
        opt.handles.hq = P.handles.hq;
        opt.handles.hs = P.handles.hs;
        opt.handles.ht = P.handles.ht;
        opt.query = 'replace';
    catch
        P.handles.axesECD = axes(...
            'parent',P.handles.hfig,...
            'Position',[0.13 0.55 0.775 0.4],...
            'hittest','off',...
            'visible','off',...
            'deleteFcn',@back2defaults);
        opt.ParentAxes = P.handles.axesECD;
        opt.hfig = P.handles.hfig;
    end
    w = reshape(mu_w,3,[]);
    s = reshape(mu_s, 3, []);
    [out] = spm_eeg_displayECD(...
        s,w,reshape(diag(S_s),3,[]),[],opt);
        P.handles.hp = out.handles.hp;
        P.handles.hq = out.handles.hq;
        P.handles.hs = out.handles.hs;
        P.handles.ht = out.handles.ht;
end

% plot data and predicted data
pos = P.forward.sens.prj;
ChanLabel = P.channels;
in.f = P.handles.hfig;
in.noButtons = 1;
try
    P.handles.axesY;
catch
    figure(P.handles.hfig)
    P.handles.axesY = axes(...
        'Position',[0.02 0.3 0.3 0.2],...
        'hittest','off');
    in.ParentAxes = P.handles.axesY;
    spm_eeg_plotScalpData(y,pos,ChanLabel,in);
    title(P.handles.axesY,'measured data')
end
if isempty(flag) || isequal(flag,'data') || isequal(flag,'ecd')
    yHat = P.gmn*mu_w;
    miY = min([yHat;y]);
    maY = max([yHat;y]);
    try
        P.handles.axesYhat;
        d = get(P.handles.axesYhat,'userdata');
        yHat = yHat(d.goodChannels);
        clim = [min(yHat(:))-( max(yHat(:))-min(yHat(:)) )/63,...
            max(yHat(:))];
        ZI = griddata(...
            d.interp.pos(1,:),d.interp.pos(2,:),full(double(yHat)),...
            d.interp.XI,d.interp.YI);
        set(d.hi,'Cdata',flipud(ZI));
        caxis(P.handles.axesYhat,clim);
        delete(d.hc)
        [C,d.hc] = contour(P.handles.axesYhat,flipud(ZI),...
            'linecolor',0.5.*ones(3,1));
        set(P.handles.axesYhat,...
            'userdata',d);
    catch
        figure(P.handles.hfig)
        P.handles.axesYhat = axes(...
            'Position',[0.37 0.3 0.3 0.2],...
            'hittest','off');
        in.ParentAxes = P.handles.axesYhat;
        spm_eeg_plotScalpData(yHat,pos,ChanLabel,in);
        title(P.handles.axesYhat,'predicted data')
    end
    try
        P.handles.axesYhatY;
    catch
        figure(P.handles.hfig)
        P.handles.axesYhatY = axes(...
            'Position',[0.72 0.3 0.25 0.2],...
            'NextPlot','replace',...
            'box','on');
    end
    plot(P.handles.axesYhatY,y,yHat,'.')
    set(P.handles.axesYhatY,...
        'nextplot','add')
    plot(P.handles.axesYhatY,[miY;maY],[miY;maY],'r')
    set(P.handles.axesYhatY,...
        'nextplot','replace')
    title(P.handles.axesYhatY,'predicted vs measured data')
    axis(P.handles.axesYhatY,'square','tight')
    grid(P.handles.axesYhatY,'on')

end

if isempty(flag) || isequal(flag,'var')
    % plot precision hyperparameters
    try
        P.handles.axesVar1;
    catch
        figure(P.handles.hfig)
        P.handles.axesVar1 = axes(...
            'Position',[0.05 0.05 0.25 0.2],...
            'NextPlot','add',...
            'box','on');
    end
    plot(P.handles.axesVar1,it,log(a1./b1),'.')
    title(P.handles.axesVar1,'measurement noise precision (log)')
    axis(P.handles.axesVar1,'square','tight')
    grid(P.handles.axesVar1,'on')

    try
        P.handles.axesVar2;
    catch
        figure(P.handles.hfig)
        P.handles.axesVar2 = axes(...
            'Position',[0.37 0.05 0.25 0.2],...
            'NextPlot','add',...
            'box','on');
    end
    plot(P.handles.axesVar2,it,log(a2./b2),'.')
    title(P.handles.axesVar2,'ECD orientations precision (log)')
    axis(P.handles.axesVar2,'square','tight')
    grid(P.handles.axesVar2,'on')

    try
        P.handles.axesVar3;
    catch
        figure(P.handles.hfig)
        P.handles.axesVar3 = axes(...
            'Position',[0.72 0.05 0.25 0.2],...
            'NextPlot','add',...
            'box','on');
    end
    plot(P.handles.axesVar3,it,log(a3./b3),'.')
    title(P.handles.axesVar3,'ECD location precision (log)')
    axis(P.handles.axesVar3,'square','tight')
    grid(P.handles.axesVar3,'on')

end

if ~isempty(flag) && (isequal(flag,'ecd') || isequal(flag,'mGN') )
    try
        P.handles.hte(2);
    catch
        figure(P.handles.hfig)
        P.handles.hte(2) = uicontrol('style','text',...
            'units','normalized',...
            'position',[0.2,0.91,0.6,0.02],...
            'backgroundcolor',[1,1,1]);
    end
    set(P.handles.hte(2),'string',...
        ['ECD locations: Modified Gauss-Newton scheme... ',num2str(floor(P.pc)),'%'])
else
    try
        set(P.handles.hte(2),'string','VB updates on hyperparameters')
    end       
end

try
    P.handles.hte(1);
catch
    figure(P.handles.hfig)
    P.handles.hte(1) = uicontrol('style','text',...
        'units','normalized',...
        'position',[0.2,0.94,0.6,0.02],...
        'backgroundcolor',[1,1,1]);
end
try
    set(P.handles.hte(1),'string',...
        ['Model evidence: p(y|m) >= ',num2str(F(end),'%10.3e\n')])
end

try
    P.handles.hti;
catch
    figure(P.handles.hfig)
    P.handles.hti = uicontrol('style','text',...
        'units','normalized',...
        'position',[0.3,0.97,0.4,0.02],...
        'backgroundcolor',[1,1,1],...
        'string',['VB ECD inversion: trial #',num2str(P.ltr(P.ii))]);
end
drawnow

function back2defaults(e1,e2)
hf = spm_figure('FindWin','Graphics');
P = get(hf,'userdata');
try
    set(hf,'colormap',P.handles.SPMdefaults.col);
    set(hf,'renderer',P.handles.SPMdefaults.renderer);
end





