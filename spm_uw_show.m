function spm_uw_show(mode,p1,p2,p3,p4,p5,p6)
% Manage graphical output for spm_uw_estimate
%
% FORMAT spm_uw_show(mode,p1,...)
%
% mode      - Verb specifying action.
% p1-p6     - Depends on mode.
%
% FORMAT spm_uw_show('FinIter',SS,beta,fot,sot,ref,q)
%
%__________________________________________________________________________
% Copyright (C) 2002-2013 Wellcome Trust Centre for Neuroimaging

% Jesper Andersson
% $Id: spm_uw_show.m 5645 2013-09-19 17:58:45Z guillaume $


persistent iter;

eff_string = {'x-trans','y-trans','z-trans','pitch','roll','yaw'};

switch mode
   case 'SmoothStart' % Report on start of smoothing.
      spm_progress_bar('Init',p1,...
      'Creating temporary smooth files',...
      'Scans completed');
   case 'SmoothUpdate'
      spm_progress_bar('Set',p1);
   case 'SmoothEnd'
      spm_progress_bar('Clear');
   case 'MaskStart' % Report on start of making mask.
      spm_progress_bar('Init',p1,...
      'Ascertaining data in all voxels and scans',...
      'Scans completed');
   case 'MaskUpdate'
      spm_progress_bar('Set',p1);
   case 'MaskEnd'
      spm_progress_bar('Clear');
   case 'NewIter'     % Report on start of new iteration.
      iter = p1;
   case 'StartRef'    % Report on start of evaluation of effects on reference scan.
      spm_progress_bar('Init',p1,...
      sprintf('Computing effects on reference scans, iteration %d',iter),...
      'Scans completed');
   case 'NewRef'      % Report on new reference scan.
      spm_progress_bar('Set',p1);
   case 'EndRef'
      spm_progress_bar('Clear');
   case 'StartAtA'    % Report on start of evaluation of AtA.
      spm_progress_bar('Init',p1,sprintf('Computing AtA, iteration %d',iter),...
      'Sub-matrices completed');
   case 'NewAtA'      % Report on new scan for AtA.
      spm_progress_bar('Set',p1);
   case 'EndAtA'
      spm_progress_bar('Clear');
   case 'StartAty'    % Report on start of evaluation of Aty.
      spm_progress_bar('Init',p1,sprintf('Computing Aty, iteration %d',iter),...
      'Scans completed');
   case 'NewAty'      % Report on new scan for Aty.
      spm_progress_bar('Set',p1);
   case 'EndAty'
      spm_progress_bar('Clear');
   case 'StartInv'    % Report on start of inversion of AtA.
      spm_progress_bar('Init',p1,sprintf('Inverting AtA, iteration %d',iter),'');
   case 'EndInv'
      spm_progress_bar('Clear');
   case 'StartRefit'
      spm_progress_bar('Init',p1,'Refitting field on regular grid','');
   case 'EndRefit'
      spm_progress_bar('Clear');
   case 'FinIter'     % Report on outcome of iteration.
              % p1 = Residual sum of squares.
              % p2 = Current derivative fields
              % p3 = fot
              % p4 = sot
              % p5 = Any scan
              % p6 = Matrix of mean corrected movement parameters.
      if nargin ~= 7
         error('Wrong no. of arguments in call to ud_graphwip');
      end
      fg = spm_figure('FindWin','Graphics');
      if ~isempty(fg)
          spm_figure('Clear','Graphics');
          %
          % Get mask for display of deformation maps.
          %
          mask = zeros(p5.dim(1:3));
          spm_smooth(reshape(p5.dat,p5.dim(1:3)),mask,5);
          mask = mask > 0.3*mean(mean(mean(mask)));
       
          %
          % Display deformation maps.
          %
          spm_figure('Select', fg);
          fg_pos = get(fg,'Position');
          fs = spm('FontSizes');
          pf = spm_platform('fonts');

          %
          % Write general title.
          %
          uicontrol(fg,'Style','Text',...
                    'Position',[fg_pos(3)/10 0.95*fg_pos(4) 8*fg_pos(3)/10 30],...
                    'String','Estimation of EPI deformation fields',...
                    'ForegroundColor','k','BackgroundColor','w',...
                    'FontSize',fs(20),'FontName',pf.times)

          %
          % Code modelled on spm_check_registration.m
          %          
          mn = size(p2,2);
          n  = round(mn^0.4);
          m  = ceil(mn/n);
          w  = 1/n;
          h  = .75/m;
          ds = (w+h)*0.02;

          spm_orthviews('Reset');
          global st;
          P = p5; 
          for ij=1:mn
             P.dat = p2(:,ij);
             P.dat = reshape(P.dat,P.dim(1:3)) .* mask;
             i  = 1-h*(floor((ij-1)/n)+1);
             j  = w*rem(ij-1,n);
             gh = spm_orthviews('Image',P,[j+ds/2 i+ds/2 w-ds h-ds]);
             if ij==1 
                 spm_orthviews('Space');
                 spm_orthviews('maxBB'); 
             end
             if ij <= size(p3,2)  % If first order effect. 
                 string = sprintf('Derivative w.r.t. %s',eff_string{p3(ij)});
             else
                 string = sprintf('Second order effect w.r.t. %s and %s',...
                     eff_string{p4(ij-size(p3,2),1)},eff_string{p4(ij-size(p3,2),2)});
             end
             if mn > 6 
                 fsi = 12;
             elseif mn > 4
                 fsi = 14;
             else
                 fsi = 16;
             end
             uicontrol(fg,'Style','Text',...
                    'Position',[(st.vols{gh}.area(1)+(st.vols{gh}.area(3)+ds)/2)*fg_pos(3)...
                                st.vols{gh}.area(2)*fg_pos(4)...
                                ((st.vols{gh}.area(3)-ds)/2)*fg_pos(3)...
                                (2*st.vols{gh}.area(4)/5)*fg_pos(4)],...
                    'String',string,...
                    'ForegroundColor','k','BackgroundColor','w',...
                    'FontSize',fs(fsi),'FontName',pf.times)
          end
          %
          % Show plot of residual squared error
          %
          ax = axes('Position',[.1 .025 .375 .17]);
          indx = find(p1 ~= 0); 
          plot(ax,indx,p1(indx),'-k','LineWidth',2);
          title(ax,'Residual error','FontSize',fs(14),'FontName',pf.times);

          %
          % Show plot of relevant movement parameters.
          %
          ax = axes('Position',[.575 .025 .375 .17]);
          plot(ax,p6);
          title(ax,'Movement parameters','FontSize',fs(14),'FontName',pf.times);
      end
             
   case 'FinTot'     % Report on outcome of undeformation.
      spm_progress_bar('Clear');
   otherwise         % Ignore unknown actions.

drawnow;

end
