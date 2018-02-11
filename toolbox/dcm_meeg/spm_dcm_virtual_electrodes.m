function CSD = spm_dcm_virtual_electrodes(DCM,s,p,TYPE)
% Posterior estimates of coupling among selected populations
% FORMAT CSD = spm_dcm_virtual_electrodes(DCM,s,p,TYPE)
%
% DCM  -  inverted DCM (see below)
% s    -  indices of source, node or region
% p    -  indices of population in node
% TYPE -  {'CSD','CCF','DTF','GCA','COH','FSD'}
%
% If called  with an output argument, graphics will be suppressed
%
% Example:
% >> spm_dcm_virtual_electrodes(DCM,[1,2,1],[1,1,8],'DTF')
%
% Estimates:
%--------------------------------------------------------------------------
% DCM.dtf                   - directed transfer functions (source space)
% DCM.ccf                   - cross covariance functions (source space)
% DCM.coh                   - cross coherence functions (source space)
% DCM.fsd                   - specific delay functions (source space)
% DCM.pst                   - peristimulus time
% DCM.Hz                    - frequency
%
% DCM.Ep                    - conditional expectation
% DCM.Cp                    - conditional covariance
% DCM.Pp                    - conditional probability
% DCM.Hc                    - conditional responses (y), channel space
% DCM.Rc                    - conditional residuals (y), channel space
% DCM.Hs                    - conditional responses (y), source space
% DCM.Ce                    - ReML error covariance
% DCM.F                     - Laplace log evidence
% DCM.ID                    - data ID
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_virtual_electrodes.m 7256 2018-02-11 14:45:00Z karl $

% population specific cross spectra
%--------------------------------------------------------------------------
M     = rmfield(DCM.M,'U');
ns    = length(DCM.Sname);        % number of sources
ne    = length(s);                % number of electrodes
qp    = DCM.Ep;                   % conditional expectations
qp.b  = zeros(2,ne) - 32;         % suppress non-specific and
qp.c  = zeros(2,ne) - 32;         % specific channel noise
M.l   = ne;

% defaults
%--------------------------------------------------------------------------
if nargin < 4, TYPE = 'CSD';             end

% specify the j-th population in the i-th source
%--------------------------------------------------------------------------
if isnumeric(qp.J)
    J      = cell(1,ns);
    [J{:}] = deal(qp.J);
    qp.J   = J;
end
try
    for i = 1:ne
        J             = spm_zeros(qp.J);
        J{s(i)}(p(i)) = 1;
        qp.LG(i,:)    = spm_vec(J)';
    end
catch
    disp('Please check you indices:')
    disp('Sources:')
    disp(DCM.Sname)
    return
end


% get Jacobian and cross spectral density
%--------------------------------------------------------------------------
[~,~,~,~,~,dfdx] = spm_dcm_mtf(qp,M);
[csd,Hz,dtf,gew] = spm_csd_mtf(qp,M,DCM.xU);  % conditional cross spectra

% set metric, support and labels/strings
%--------------------------------------------------------------------------
nc   = length(csd);                  % number of conditions or trials
switch TYPE
    
    case('CSD') % conditional cross spectral density
        %------------------------------------------------------------------
        for i = nc, CSD{i} = abs(csd{i}); end
        X    = Hz;
        xstr = 'Frequency (Hz)';
        ystr = 'Cross spectra';
        
    case('CCF') % conditional covariance functions
        %------------------------------------------------------------------
        [CSD,X] = spm_csd2ccf(csd,Hz);
        xstr = 'Lag (sec)';
        ystr = 'Covariance';
        
    case('DTF') % conditional covariance functions
        %------------------------------------------------------------------
        for i = nc, CSD{i} = abs(dtf{i}); end
        X    = Hz;
        xstr = 'Frequency (Hz)';
        ystr = 'DTF';
        
    case('GCA') % conditional covariance functions
        %------------------------------------------------------------------
        CSD  = gew;
        X    = Hz;
        xstr = 'Frequency (Hz)';
        ystr = 'Granger';
        
    case('COH') % conditional covariance functions
        %------------------------------------------------------------------
        CSD  = spm_csd2coh(csd,Hz);
        X    = Hz;
        xstr = 'Frequency (Hz)';
        ystr = 'Coherence';

    case('FSD') % conditional covariance functions
        %------------------------------------------------------------------
        [coh,CSD] = spm_csd2coh(csd,Hz);
        X    = Hz;
        xstr = 'Frequency (Hz)';
        ystr = 'Delay';
        
    otherwise
        
        disp('please specify TYPE of coupling')
        disp({'CSD','CCF','DTF','GCA','COH','FSD'})
        return
        
end

% return
%--------------------------------------------------------------------------
if nargout, return, end

% graphics
%--------------------------------------------------------------------------
spm_figure('GetWin','Virtual electrodes'); clf
co   = {'b', 'r', 'g', 'm', 'y', 'k', 'c'};
name = DCM.Sname(s);

% for each pair of populations
%--------------------------------------------------------------------------
for i = 1:ne
    for j = 1:ne
        subplot(ne + 1,ne,(i - 1)*ne + j)
        for k = 1:nc
            str = sprintf('%s: %s to %s',ystr,name{j},name{i});
            plot(X,CSD{k}(:,i,j),'color',co{k}), hold on
            title(str,'FontSize',14,'FontWeight','Bold')
            xlabel(xstr),ylabel(ystr)
            axis square, spm_axis tight
        end
    end
end
str   = cell(1,nc);
for k = 1:nc, str{k} = sprintf('trial: %i',k); end
legend(str)

% Effective connectivity (Jacobian)
%--------------------------------------------------------------------------
subplot(ne + 1,1,ne + 1)
imagesc(dfdx), axis square
title('Jacobian','FontSize',14,'FontWeight','Bold')
xlabel('Hidden state'),ylabel('Hidden state')
hold on
for i = 1:ne
    j = find(qp.LG(i,:));
    plot(j,j,'o','MarkerSize',16,'color',co{i})
end
hold off
legend(name)

