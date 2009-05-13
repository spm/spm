function dipoles  = spm_eeg_dipoles_ui
% Get dipole locations and orientations from the user
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak, Christophe Phillips
% $Id: spm_eeg_dipoles_ui.m 3120 2009-05-13 13:01:03Z vladimir $

SVNrev = '$Rev: 3120 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FnUIsetup','Specify dipole model',0);

label  = {};
pnt    = [];
ori    = [];

dip_q = 0; % number of dipole 'elements' added (single or pair)
dip_c = 0; % total number of dipoles in the model
adding_dips = 1;
while adding_dips
    if dip_q>0,
        msg_dip =['Add dipoles to ',num2str(dip_c),' or stop?'];
        dip_ch = 'Single|Pair|Stop';
        dip_val = [1,2,0];
        def_opt=3;
    else
        msg_dip =['Add dipoles to model'];
        def_opt=1;
        dip_ch = 'Single|Pair';
        dip_val = [1,2];
    end
    a_dip = spm_input(msg_dip,1+dip_q,'b',dip_ch,dip_val,def_opt);
    if a_dip == 0
        adding_dips = 0;
    else
        clabel = spm_input('Source label', '+1', 's');

        if a_dip == 1
            % add a single dipole to the model
            dip_q = dip_q+1;
            % informative location prior
            str = 'Location';
     
            loc = spm_input(str, 1+dip_q+1,'e',[0 0 0]);
     
            % Moment prior
            wpr_q = spm_input('Oriented dipole?',1+dip_q+1,'b', ...
                'Yes|No',[1,0],2);
            if wpr_q
                % informative moment prior
                m = spm_input('Orientation', ...
                    1+dip_q+1,'e',[0 0 0]);
                m = m/norm(m);

                pnt    = [pnt; loc];
                ori = [ori; m];
                label       = [label; {clabel}];
            else
                pnt    =   [pnt; loc; loc; loc];
                ori =   [ori; 1 0 0; 0 1 0; 0 0 1];
                label       =   [label; {[clabel '_X']; [clabel '_Y']; [clabel '_Z']}];
            end
            dip_c = dip_c+1;
        else
            % add a pair of symmetric dipoles to the model
            dip_q = dip_q+1;
            % informative location prior
            str = 'Location (right only)';
          
            loc = spm_input(str, 1+dip_q+1,'e',[0 0 0]);
          
            % Moment prior
            wpr_q = spm_input('Oriented dipoles?',1+dip_q+1,'b', ...
                'Yes|No',[1,0],2);

            if wpr_q
                % informative moment prior
                m = spm_input('Orientation (right only)', ...
                    1+dip_q+1,'e',[0 0 0]);
                m = m/norm(m);

                pnt    = [pnt; loc; -loc(1) loc(2:3)];
                ori = [ori; m; -m(1) m(2:3)];
                label       = [label; {[clabel 'R']; [clabel 'L']}];
            else
                pnt    =   [pnt; repmat(loc, 3, 1); repmat([-loc(1) loc(2:3)], 3, 1)];
                ori =   [ori; 1 0 0; 0 1 0; 0 0 1; 1 0 0; 0 1 0; 0 0 1];
                label       =   [label; {[clabel 'R_X']; [clabel 'R_Y']; [clabel 'R_Z']; [clabel 'L_X']; [clabel 'L_Y']; [clabel 'L_Z']}];
            end

            dip_c = dip_c+2;
        end
    end
end

spm_figure('GetWin','Graphics'); clf

sdip= [];
sdip.n_seeds = 1;
sdip.n_dip   = size(pnt, 1);
sdip.Mtb     = 1;
sdip.j{1}    = reshape(ori', [], 1);
sdip.loc{1}  = pnt';
spm_eeg_inv_ecd_DrawDip('Init', sdip)
spm_eeg_inv_ecd_DrawDip('drawdip', 1, size(pnt, 1)+1);


dipoles = [];
dipoles.pnt   = pnt;
dipoles.ori   = ori;
dipoles.label = label;
