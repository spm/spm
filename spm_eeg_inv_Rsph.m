function [spheres,dipoles,L,Lan,flags] = spm_eeg_inv_Rsph(model,dipoles,flags)

%____________________________________________________________________________
%
% FUNCTION [spheres,dipoles,L,Lan,flags] = spm_eeg_inv_Rsph(model,dipoles,flags)
% or [spheres] = spm_eeg_inv_Rsph(model,[])
%
% Estimate the best spherical model for the realistic head (scalp or brain) model.
% Then, if a dipoles structure is provided, it deforms the dipoles
% distribution to fit the sphere model and calculate the leadfield
% for the model & dipoles, if required.
%
% The model could only contain 1 volume model, brain or scalp.
% It assumes the volume is the scalp, if 'br_only' is not a field of 'model' or
% if 'model.br_only' equals 0.
% When only the brain is used, the skull and scalp thickness can be provided through
% the field 'sk_sc_th' in model (skull thickness first)
%
% IN :
%   - model : head model with the electrode setup. It must contain the mesh
%               for the scalp at the least.
%   - dipoles : dipoles structure, could be cortical mesh or regular cubic grid
%   - flags   : useful flags for the routine
%       * br_only   : model contins only the brain volume, as if based on PET/EPI
%       * use_br    : use the brain surface to define the sphere model (1) or not (0)
%       * figs      : display colourful images of the model (1) or not (0)
%       * sk_sc_th  : a priori skull and scalp thickness [8 8]mm
%       * calc_dip  : deal with dipoles structure (1) or not (0)
%       * calc_L    : calculate the leadfield for model (with dipoles) (1) or not (0)
% OUT :
%   - spheres : new structure 'sphere' for spherical model.
%   - dipoles : dipoles structure with new structure 'dipS' for "sphere dipoles".
%   - L       : leadfiled for dipoles, oriented if the orientation is provided.
%   - Lan     : leadfiled for dipoles, orientation free.
%   - flags   : flags used.
%
% model : same structure as input but the extra sphere structure
% .spheres
%   .centre   : centre of the sphere model (in mm)
%   .Rsc \
%   .Rsk  |   : radii of the scalp, skull and brain sphere
%   .Rbr /
%   .r_s      : original radius of the un-deformed "spheres" (scalp/brain surface)
%   .Ref_surf : Scalp(1) or brain(2), default 1
%   .R        : Radius of reference surface
%
% dipS : same structure as dipoles (no L or R part though) with some extra stuff
%   .centre : centre of the sphere model (in mm)
%   .Rsc \
%   .Rsk  |   : radii of the scalp, skull and brain sphere
%   .Rbr /
%   .r_s      : original radius of the un-deformed "spheres" (scalp/brain surface)
%   .Ref_surf : Scalp(1) or brain(2), default 1
%   .r_d    : original radius of the dipoles set
%   .r_ds   : 'radius' of scalp in direction of each dipole
%
% The spherical approximation comes from the following paper:
% Laurent Spinelli, Sara Gonzalez Andino, Goran Lantz, Margitta Seeck, Christoph M. Michel
% Electromagntic inverse solutions in anatomically constrained spherical head models
% Brain Topography, vol. 13, 2 2000
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Christophe Phillips,
% $Id: spm_eeg_inv_Rsph.m 1039 2007-12-21 20:20:38Z karl $

def_flags = struct('br_only',0,'use_br',0,'figs',0,'sk_sc_th',[8 8],...
    'calc_dip',0,'calc_L',0);

if nargin<3
    flags = def_flags;
else
    fnms  = fieldnames(def_flags);
    for i = 1:length(fnms),
        if ~isfield(flags,fnms{i}),
            flags = setfield(flags,fnms{i},getfield(def_flags,fnms{i}));
        end
    end
end

% Pd : filename of dipole model file
% Pm : filename of head model file

if nargin  < 2
    Pd = spm_select(Inf,'*dipoles*.mat','Dipole set or nothing');
    if isempty(Pd)
        flags.calc_dip = 0;
        flags.calc_L = 0;
    else
        load(Pd)
        flags.calc_dip = 1;
    end
else
    Pd = [];
    if ~isempty(dipoles)
        flags.calc_dip = 1;
        flags.calc_L = 1;
    end
end
Pm = [];
if nargin < 1
    Pm = spm_select(1,'*model*.mat','Head model');
    load(Pm)
end

if ~isfield(model,'spheres')
    
    % Deal with the head model first, if it isn't already done !
    %===============================
    % Scalp & brain vertices in mm + Electrodes location in mm
    fprintf('\n\n'), fprintf('%c','='*ones(1,80)), fprintf('\n')
    fprintf(['\tGenerate realistic sphere model from tesselated surfaces.\n']);

    Nsurf = length(model.head) ;
    if Nsurf==3 | Nsurf==2  % multiple surfaces
        sXYZ = model.head(end).XYZmm ;
        bXYZ = model.head(1).XYZmm ;
        if flags.use_br
            rXYZ = bXYZ;
            rtri = model.head(1).tri;
        else
            rXYZ = sXYZ;
            rtri = model.head(end).tri;
        end
    else    % only 1 surface, scalp or brain ?
        if isfield(model,'br_only')
            flags.br_only = model.br_only;
            if flags.br_only
                bXYZ = model.head(1).XYZmm ;
                rXYZ = bXYZ;
                flags.use_br = 1;
            else
                sXYZ = model.head(1).XYZmm ;
                rXYZ = sXYZ;
            end
        else
            sXYZ = model.head(1).XYZmm ;
            rXYZ = sXYZ;
        end
        rtri = model.head(1).tri;
    end
%     if flags.figs,
%         % Display the scalp/brain surface
%         fvb.vertices = rXYZ'; fvb.faces = rtri';
%         spm_eeg_inv_displTes(fvb)
%         if flags.use_br
%             title('Brain surface')
%         else
%             title('Scalp surface')
%         end
%     end

    if flags.use_br
        Ref_surf = 2;
    else
        Ref_surf = 1;
    end
    eXYZ = model.electrodes.XYZmm ;

    % Define the orientation of the xy plane, defining the upper part of the head
    [yM,iM] = max(eXYZ(2,:));
    [zm,im] = min(eXYZ(3,:));
    ym   = eXYZ(2,im);
    c_th = (eXYZ(3,im)-eXYZ(3,iM))/norm(eXYZ(1:2,im)-eXYZ(1:2,iM));
    s_th = sin(acos(c_th));
    no = [c_th s_th] ;
    no = no/norm(no);
    if no(2)<0, no = -no; end % Make sure no points upward
    d = 20;


    l_use  = find((no*rXYZ(2:3,:))>=-d);
    l_Nuse = 1:size(rXYZ,2);
    l_Nuse(l_use) = [];

    if flags.figs & flags.use_br
        spm_eeg_inv_displScEl(model.head(1),model.electrodes);
    else
        spm_eeg_inv_displScEl(model.head(end),model.electrodes);
    end
    hold on, plot3(rXYZ(1,l_Nuse),rXYZ(2,l_Nuse),rXYZ(3,l_Nuse),'bo')
    rotate3d on
    title('Ref. surface with electrodes, blue vertices are not considered')

    P      = [rXYZ(:,l_use)' ones(length(l_use),1)]\sum(rXYZ(:,l_use).^2)' ;
    centre = P(1:3)/2;
    R      = sqrt(P(4)+sum(centre.^2));
    if flags.use_br
        fprintf('Brain sphere will have centre : %3.2f %3.2f %3.2f, and radius R= %3.2f \n', ...
            [centre' R])
    else
        fprintf('Scalp sphere will have centre : %3.2f %3.2f %3.2f, and radius R= %3.2f \n', ...
            [centre' R])
    end

    % Sphere defintion for brain: same centre but Rbr
    if Nsurf > 1
        if flags.use_br
            l_use = find((no*sXYZ(2:3,:))>=0);
            Rbr = R;
            Rsc = sqrt(mean(sum((sXYZ(1:3,l_use)-centre*ones(1,length(l_use))).^2)));
            Rsk = (R+Rsc)/2;
        else
            l_use = find((no*bXYZ(2:3,:))>=0);
            Rbr = sqrt(mean(sum((bXYZ(1:3,l_use)-centre*ones(1,length(l_use))).^2)));
            Rsk = (R+Rbr)/2;
            Rsc = R;
        end
    else
        if flags.use_br
            Rbr = R;
            Rsk = Rbr+flags.sk_sc_th(1);
            Rsc = Rsk+flags.sk_sc_th(2);
        else
            Rsc = R ;
            Rsk = Rsc-flags.sk_sc_th(2);
            Rbr = Rsk-flags.sk_sc_th(1);
        end
    end
    fprintf('Radii for scalp: %3.2f, skull: %3.2f & brain: %3.2f spheres.\n', ...
            [Rsc Rsk Rbr])

    % Scalp(brain) vertices in spherical coordinates
    % (azimuth TH, elevation PHI, and radius R)
    [th,phi,r_s] = cart2sph(rXYZ(1,:)-centre(1),rXYZ(2,:)-centre(2), ...
                            rXYZ(3,:)-centre(3));
    % Scalp/brain vertices *on* a sphere of radius 1, in Cartesian coord
    [Xss,Yss,Zss] = sph2cart(th, phi, 1);
    SsXYZ = [Xss;Yss;Zss] ;

    if flags.figs
        % Display "spherical" scalp, with colour coded radius.
        fv_sp.vertices = SsXYZ'; fv_sp.faces = rtri';
        spm_eeg_inv_displTes(fv_sp,r_s',1)
        title('"Spherical" scalp/brain, with colour coded radius.')
    end

    if flags.use_br
        % Avoid that the lower part of the brain gets too stretched
        Mr_use = 1.5*max(abs(R-r_s(l_use)));
        l_too_short = find(r_s(l_Nuse)<(R-Mr_use));
        r_s(l_Nuse(l_too_short)) = R-Mr_use;
    end

    % Electrodes in spherical coordinates (azimuth TH, elevation PHI, and radius R)
    [th,phi,r_el] = cart2sph(eXYZ(1,:)-centre(1),eXYZ(2,:)-centre(2),eXYZ(3,:)-centre(3));
    % Electrodes *on* the scalp sphere in Cartesian coord
    [Xes,Yes,Zes] = sph2cart(th, phi, Rsc);
    SseXYZ = [Xes;Yes;Zes] ;

elseif flags.calc_dip
    centre   = model.spheres.centre;
    SsXYZ    = model.spheres.Ss_XYZ;
    r_s      = model.spheres.r_s;
    R        = model.spheres.R;
    Ref_surf = model.spheres.Ref_surf;
else
    disp('I''ve nothing to do here. Thank you for calling...')
    return
end

if flags.calc_dip
    % Then take care of the dipoles
    %==============================
    % Decide if I'm dealing with a cortical mesh or a regular grid
    if isfield(dipoles,'vert')
        cort_or_grid = 1;   % Cortical mesh
    elseif isfield(dipoles,'XYZmm')
        cort_or_grid = 2;   % 3D grid
    else
        error('Don''t know this type of dipoles structure')
    end
    if cort_or_grid==1
        [th,phi,r_d] = cart2sph(dipoles.vert(1,:)-centre(1), ...
            dipoles.vert(2,:)-centre(2), ...
            dipoles.vert(3,:)-centre(3));
    else
        XYZmm = dipoles.XYZmm;
        [th,phi,r_d] = cart2sph(XYZmm(1,:)-centre(1), ...
            XYZmm(2,:)-centre(2), ...
            XYZmm(3,:)-centre(3));
    end
    dTPR  = [th ; phi ; r_d];
    [XSd,YSd,ZSd] = sph2cart(dTPR(1,:), dTPR(2,:), 1);
    SdXYZ = [XSd;YSd;ZSd] ;

    % Interpolation of the scalp radius for each dipoles,
    % throug a k-nearest neighbour approach.
    Nneighb = 4;
    m = -2;
    % 'radius' of scalp in direction of each dipole
    r_ds = zeros(dipoles.nr(1),1);
    for ii=1:dipoles.nr(1)
        [cos_angl,perm] = sort(abs(acos(SdXYZ(:,ii)'*SsXYZ)));
        sc_i  = perm(1:Nneighb);
        d_neighb = sqrt(sum((SdXYZ(:,ii)*ones(1,Nneighb)-SsXYZ(:,sc_i)).^2));
        r_ds(ii) = sum(d_neighb.^m .* r_s(sc_i))/sum(d_neighb.^m);
    end

    if cort_or_grid==1 & figs
        % Display cort surface with colour coded scalp radius
        fv_c.vertices = dipoles.vert'; fv_c.faces = dipoles.tri'; spm_eeg_inv_displTes(fv_c,r_ds)
    end

    % "New" dipoles coordinates after the spherical deformation
    r_dsS = r_d' .* (R./r_ds);
    [Xsd,Ysd,Zsd] = sph2cart(dTPR(1,:), dTPR(2,:), r_dsS');
    SdXYZ = [Xsd;Ysd;Zsd];

    % Create a new structure for the dipoles in "spherical" head.
    if cort_or_grid==1
        dipS = struct('vert',[SdXYZ],'tri',dipoles.tri,'or',[], ...
            'Nv',dipoles.Nv,'Nf',dipoles.Nf, ...
            'nr',[dipoles.Nv dipoles.Nf],'centre',centre, ...
            'Ref_surf',Ref_surf,...
            'R',R,'Rsc',[],'Rsk',[],'Rbr',[],'r_d',r_d','r_ds',r_ds);
    else
        dipS = struct('XYZmm',[SdXYZ],'or',[],'nr',dipoles.nr,'centre',centre, ...
            'Ref_surf',Ref_surf,...
            'R',R,'Rsc',[],'Rsk',[],'Rb',[],'r_d',r_d','r_ds',r_ds);
    end


    if cort_or_grid==1 & figs
        % Display deformed cort surface
        fv_Sc.vertices = dipS.vert'; fv_Sc.faces = dipoles.tri';
        spm_eeg_inv_displTes(fv_Sc,(R./r_ds),1);
    end

    if cort_or_grid==1
        % Estimate the orientation in the new configuration
        t_or = zeros(3,dipoles.Nv);
        unit = [1 1 1];
        for ii=1:dipoles.Nf
            va = dipoles.vert(:,dipoles.tri(2,ii)) - ...
                dipoles.vert(:,dipoles.tri(1,ii));
            vb = dipoles.vert(:,dipoles.tri(3,ii)) - ...
                dipoles.vert(:,dipoles.tri(1,ii));
            vn = [va(2)*vb(3)-va(3)*vb(2) va(3)*vb(1)-va(1)*vb(3) va(1)*vb(2)-va(2)*vb(1)];
            vn = vn'/norm(vn) ;
            t_or(:,dipoles.tri(:,ii)) = t_or(:,dipoles.tri(:,ii))+vn*unit;
        end
        dipS.or = t_or./( unit'*sqrt(sum((t_or).^2)) );
    else
        if isfield(dipoles,'or')
            dipS.or = dipoles.or;
        else
            dipS.or = [];
        end
    end
end

% Check the radius of the brain sphere if dipoles were provided
%==============================================================
if flags.calc_dip
    %   Make sure all dipoles are in the brain sphere
    Rb_ch = max(Rbr,ceil(max(r_dsS+2)));
else
    Rb_ch = Rbr;
end

if Rbr~=Rb_ch
    fprintf('Problem with dipoles outside the brain volume\n')
    if flags.br_only
        Rsk_ch = Rb_ch+model.sk_sc_th(1);
        Rsc_ch = Rsk_ch+model.sk_sc_th(2);
        fprintf('\tChanging scalp radius from %3.1g to %3.1g,\n',Rsc,Rsc_ch)
        fprintf('\t     and skull radius from %3.1g to %3.1g,\n',Rsk,Rsk_ch)
        Rbr = Rb_ch ; Rsk = Rsk_ch ; Rsc = Rsc_ch ;
    else
        Rsk_ch = (Rsc+Rb_ch)/2;
        fprintf('\tChanging brain radius from %3.2f to %3.2f,\n',Rbr,Rb_ch)
        fprintf('\t     and skull radius from %3.2f to %3.2f,\n',Rsk,Rsk_ch)
        Rbr = Rb_ch ; Rsk = Rsk_ch ;

    end
end
if flags.calc_dip
    dipS.Rbr = Rbr;
    dipS.Rsc = Rsc;
    dipS.Rsk = Rsk;
    fprintf('Skull radius Rsk= %3.2f and brain radius Rb= %3.2f.\n',[Rsk Rbr])
    dipoles.dipS = dipS;
else
    dipoles = 'not processed here';
end

% Estimate the Leadfield, if required
%=======================
if flags.calc_L
    if ~isfield(model,'param') || ~isfield(model.param,'sigma')
        model.param.sigma = [.33 .004 .33];
        fprintf('\nWarning : conductivity values not defined in the model !\n')
        fprintf(' I set them my self to %1.3g %1.3g %1.3g .\n',model.param.sigma)
    end
    % Estimate the leadfield for the "spherical" model & dipoles
    %============================================================
    fprintf('Start Leadfield calculation \n')
    [Lan,nit] = spm_eeg_inv_Lana(SdXYZ,SseXYZ,Rsc,Rsk,Rbr,model.param.sigma);
    if ~isempty(dipS.or)
        fprintf('Calculating Leadfiled of oriented dipoles.\n')
        % Oriented dipoles...
        Or1 = kron(spdiags(dipS.or(1,:)',0,dipS.nr(1),dipS.nr(1)), [1 0 0]') ;
        Or2 = kron(spdiags(dipS.or(2,:)',0,dipS.nr(1),dipS.nr(1)), [0 1 0]') ;
        Or3 = kron(spdiags(dipS.or(3,:)',0,dipS.nr(1),dipS.nr(1)), [0 0 1]') ;
        Or = Or1+Or2+Or3 ;
        L = Lan*Or;
%       Lorig = L;
    else
        fprintf('No orientation provided for the dipoles.\n')
        L = Lan;
    end
    L = L - ones(size(L,1),1)*mean(L);
else
    L   = 'Not calculated yet.';
    Lan = L;
end

if ~isempty(Pm)
    Pm_s = [spm_str_manip(Pm,'r'),'_3sph.mat'];
    fprintf('\n Saving model into file : %s\n',Pm_s)
    save(Pm_s,'model')
end
if ~isempty(Pd) %usePd
    Pd_s = [spm_str_manip(Pd,'r'),'_3sph.mat'];
    fprintf('\n Saving dipoles into file : %s\n',Pd_s)
    save(Pm_s,'dipoles')
end

spheres = struct('centre',centre,'R',R,'Rsc',Rsc,'Rsk',Rsk,'Rbr',Rbr,'r_s',r_s, ...
    'Sc_elXYZ',SseXYZ,'Ss_XYZ',SsXYZ,'Ref_surf',Ref_surf,'flags',flags);

return
