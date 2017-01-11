function [Y,y,beta,Bcov,G] = spm_graph(SPM,XYZ,xG)
% Return adjusted data for a given voxel location
% FORMAT [Y,y,beta,Bcov,G] = spm_graph(SPM,XYZ,xG)
%
% SPM    - structure containing generic details about the analysis
% XYZ    - [x y z]' coordinates {voxel}
% xG     - structure containing details about action to perform
%   .def - string describing data type to be returned. One of:
%           'Contrast estimates and 90% C.I.'
%           'Fitted responses'
%           'Event-related responses'
%           'Parametric responses'
%           'Volterra Kernels'
%   .spec - structure containing specific details about returned data
%
% Y      - fitted   data for the selected voxel
% y      - adjusted data for the selected voxel
% beta   - parameter estimates (ML or MAP)
% Bcov   - covariance of parameter estimates (ML or conditional)
% G      - structure containing further data depending on xG details
%
% See spm_graph_ui for details.
%__________________________________________________________________________
% Copyright (C) 1996-2016 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_graph.m 6985 2017-01-11 11:51:49Z guillaume $


if nargin == 3 && isstruct(SPM) && isstruct(XYZ) && ishandle(xG)
    warning('Syntax of spm_graph changed: call spm_graph_ui instead.');
    [Y,y,beta,Bcov] = spm_graph_ui(SPM,XYZ,xG); G = [];
    return
end
if nargin < 3, [xG,G] = deal(struct([])); end

%==========================================================================
%-Extract filtered and whitened data from files
%==========================================================================
Y     = [];
y     = [];

if ismember(xG.def,{'Fitted responses','Event-related responses'})
    try
        y = spm_data_read(SPM.xY.VY,'xyz',XYZ);
    catch
        try
            % remap files in SPM.xY.P if SPM.xY.VY is no longer valid
            %--------------------------------------------------------------
            SPM.xY.VY = spm_data_hdr_read(SPM.xY.P);
            y = spm_data_read(SPM.xY.VY,'xyz',XYZ);
            
        catch
            % data has been moved or renamed
            %--------------------------------------------------------------
            choice = questdlg({'Original data have been moved or renamed',...
                'How to proceed next?'},...
                [mfilename ': data files missing...'],...
                'Specify','Search','Ignore','Ignore');
            
            switch choice
                case 'Specify'
                    [SPM.xY.P,sts] = ...
                        spm_select(numel(SPM.xY.VY),'image','Select images');
                    if ~sts
                        [Y,y,beta,Bcov] = deal([]);
                        spm('Pointer','Arrow');
                        return;
                    end
                    SPM.xY.VY = spm_data_hdr_read(SPM.xY.P);
                    for i = 1:numel(SPM.xY.VY)
                        SPM.xY.VY(i).pinfo(1:2,:) = ...
                            SPM.xY.VY(i).pinfo(1:2,:)*SPM.xGX.gSF(i);
                    end
                    y = spm_data_read(SPM.xY.VY,'xyz',XYZ);
                case 'Search'
                    SPM.xY.VY = spm_check_filename(SPM.xY.VY);
                    y = spm_data_read(SPM.xY.VY,'xyz',XYZ);
                otherwise
                    y = [];
            end
        end
    end
end

if ~isempty(y), y = spm_filter(SPM.xX.K,SPM.xX.W*y); end

%-Compute residuals
%--------------------------------------------------------------------------
if isempty(y)

    % make R = NaN so it will not be plotted
    %----------------------------------------------------------------------
    R   = NaN(size(SPM.xX.X,1),1);

else
    
    % residuals (non-whitened)
    %----------------------------------------------------------------------
    R   = spm_sp('r',SPM.xX.xKXs,y);

end


%==========================================================================
%-Get parameter and hyperparameter estimates
%==========================================================================
if ~isfield(SPM,'VCbeta') % xSPM.STAT ~= 'P'

    %-Parameter estimates:   beta = xX.pKX*xX.K*y;
    %-Residual mean square: ResMS = sum(R.^2)/xX.trRV
    %----------------------------------------------------------------------
    beta  = spm_data_read(SPM.Vbeta,'xyz',XYZ);
    ResMS = spm_data_read(SPM.VResMS,'xyz',XYZ);
    Bcov  = ResMS*SPM.xX.Bcov;

else
    % or conditional estimates with
    % Cov(b|y) through Taylor approximation
    %----------------------------------------------------------------------
    beta  = spm_data_read(SPM.VCbeta, 'xyz', XYZ);

    if isfield(SPM.PPM,'VB')
        % Get approximate posterior covariance at ic
        % using Taylor-series approximation

        % Get posterior SD beta's
        Nk = size(SPM.xX.X,2);
        for k=1:Nk
            sd_beta(k,:) = spm_data_read(SPM.VPsd(k),'xyz',XYZ);
        end

        % Get AR coefficients
        nsess = length(SPM.Sess);
        for ss=1:nsess
            for p=1:SPM.PPM.AR_P
                Sess(ss).a(p,:) = spm_data_read(SPM.PPM.Sess(ss).VAR(p),'xyz',XYZ);
            end
            % Get noise SD
            Sess(ss).lambda = spm_data_read(SPM.PPM.Sess(ss).VHp,'xyz',XYZ);
        end

        % Which block are we in ?
        % this needs updating s.t xSPM contains labels of selected voxels
        v = find((SPM.xVol.XYZ(1,:)==XYZ(1))&(SPM.xVol.XYZ(2,:)==XYZ(2))&(SPM.xVol.XYZ(3,:)==XYZ(3)));
        block_index = SPM.xVol.labels(v);
        Bcov = zeros(Nk,Nk);
        for ss=1:nsess
            % Reconstuct approximation to voxel wise correlation matrix
            post_R = SPM.PPM.Sess(ss).block(block_index).mean.R;
            if SPM.PPM.AR_P > 0
                dh = Sess(ss).a(:,1)'-SPM.PPM.Sess(ss).block(block_index).mean.a;
            else
                dh = [];
            end
            dh = [dh Sess(ss).lambda(1)-SPM.PPM.Sess(ss).block(block_index).mean.lambda];
            for i=1:length(dh)
                post_R = post_R + SPM.PPM.Sess(ss).block(block_index).mean.dR(:,:,i)*dh(i);
            end
            % Get indexes of regressors specific to this session
            scol = SPM.Sess(ss).col;
            mean_col_index = SPM.Sess(nsess).col(end)+ss;
            scol = [scol mean_col_index];

            % Reconstuct approximation to voxel wise covariance matrix
            Bcov(scol,scol) = Bcov(scol,scol) + (sd_beta(scol,1)*sd_beta(scol,1)').*post_R;
        end

    else
        Bcov     = SPM.PPM.Cby;
        for j = 1:length(SPM.PPM.l)
            l    = spm_data_read(SPM.VHp(j),'xyz',XYZ);
            Bcov = Bcov + SPM.PPM.dC{j}*(l - SPM.PPM.l(j));
        end
    end
end

%-Return if plot hasn't been defined
%--------------------------------------------------------------------------
if isempty(xG) || ~isfield(xG,'def') || isempty(xG.def)
    return;
end

%==========================================================================
%-Compute estimates
%==========================================================================
CI    = 1.6449;  % = spm_invNcdf(1 - 0.05);

switch xG.def
    
    %-Parameter estimates
    %======================================================================
    case 'Contrast estimates and 90% C.I.'
        
        Ic    = xG.spec.Ic; % contrast index
        if numel(Ic) == 1
            c = SPM.xCon(Ic).c; % contrast weights
        else
            c = Ic';
        end
        
        % compute contrast of parameter estimates and 90% C.I.
        %------------------------------------------------------------------
        cbeta = c'*beta;
        SE    = sqrt(diag(c'*Bcov*c));
        CI    = CI*SE;
        
        % returned values
        %------------------------------------------------------------------
        G.contrast      = cbeta;
        G.standarderror = SE;
        G.interval      = 2*CI;
        
        
    %-All fitted effects or selected effects
    %======================================================================
    case 'Fitted responses'
        
        Ic    = xG.spec.Ic; % contrast index
        
        % predicted or adjusted response
        %------------------------------------------------------------------
        if xG.spec.predicted

            % fitted (predicted) data (Y = X1*beta)
            %--------------------------------------------------------------
            % this should be SPM.xX.xKXs.X instead of SPM.xX.X below
            Y = SPM.xX.X*SPM.xCon(Ic).c*pinv(SPM.xCon(Ic).c)*beta;
        else

            % fitted (corrected)  data (Y = X1o*beta)
            %--------------------------------------------------------------
            Y = spm_FcUtil('Yc',SPM.xCon(Ic),SPM.xX.xKXs,beta);

        end

        % adjusted data
        %------------------------------------------------------------------
        y     = Y + R;

        % ordinate
        %------------------------------------------------------------------
        switch char(fieldnames(xG.spec.x))
            
            case 'i' % an explanatory variable

                i    = xG.spec.x.i;
                x    = SPM.xX.xKXs.X(:,i);

            case 'scan' % scan or time

                if isfield(SPM.xY,'RT')
                    x    = SPM.xY.RT*[1:size(Y,1)]';
                else
                    x    = [1:size(Y,1)]';
                end

            case 'x' % user specified

                x    = xG.spec.x.x;

        end
        
        % returned values
        %------------------------------------------------------------------
        G.x    = x;
 
        
    %-Modeling evoked responses based on Sess
    %======================================================================
    case 'Event-related responses'

        dt    = SPM.xBF.dt;
        s     = xG.spec.Sess;
        u     = xG.spec.u;
        
        % event-related response
        %------------------------------------------------------------------
        if isempty(y)
            warning(['Data not available. ' ...
                'Plotting fitted response and 90% C.I. instead.']);
            xG.spec.Rplot = 'fitted response and 90% C.I.';
        end
        switch xG.spec.Rplot
            case 'fitted response and PSTH'
                % build a simple FIR model subpartition (X); bin size = TR
                %----------------------------------------------------------
                BIN         = SPM.xY.RT;
                %BIN         = max(2,BIN);
                xBF         = SPM.xBF;
                U           = SPM.Sess(s).U(u);
                U.u         = U.u(:,1);
                xBF.name    = 'Finite Impulse Response';
                xBF.order   = round(32/BIN);
                xBF.length  = xBF.order*BIN;
                xBF         = spm_get_bf(xBF);
                BIN         = xBF.length/xBF.order;
                X           = spm_Volterra(U,xBF.bf,1);
                k           = SPM.nscan(s);
                X           = X([0:(k - 1)]*SPM.xBF.T + SPM.xBF.T0 + 32,:);

                % place X in SPM.xX.X
                %----------------------------------------------------------
                jX          = SPM.Sess(s).row;
                iX          = SPM.Sess(s).col(SPM.Sess(s).Fc(u).i);
                iX0         = [1:size(SPM.xX.X,2)];
                iX0(iX)     = [];
                X           = [X SPM.xX.X(jX,iX0)];
                X           = SPM.xX.W(jX,jX)*X;
                X           = [X SPM.xX.K(s).X0];

                % Re-estimate to get PSTH and CI
                %----------------------------------------------------------
                j           = xBF.order;
                xX          = spm_sp('Set',X);
                pX          = spm_sp('x-',xX);
                PSTH        = pX*y(jX);
                res         = spm_sp('r',xX,y(jX));
                df          = size(X,1) - size(X,2);
                bcov        = pX*pX'*sum(res.^2)/df;
                PSTH        = PSTH(1:j)/dt;
                PST         = [1:j]*BIN - BIN/2;
                PCI         = CI*sqrt(diag(bcov(1:j,(1:j))))/dt;
        end

        % basis functions and parameters
        %------------------------------------------------------------------
        X     = SPM.xBF.bf/dt;
        x     = ([1:size(X,1)] - 1)*dt;
        j     = SPM.Sess(s).col(SPM.Sess(s).Fc(u).i(1:size(X,2)));
        B     = beta(j);

        % fitted responses with standard error
        %------------------------------------------------------------------
        Y     = X*B;
        CI    = CI*sqrt(diag(X*Bcov(j,j)*X'));

        % peristimulus times and adjusted data (y = Y + R)
        %------------------------------------------------------------------
        pst   = SPM.Sess(s).U(u).pst;
        bin   = round(pst/dt);
        q     = find((bin >= 0) & (bin < size(X,1)));
        y     = R(SPM.Sess(s).row(:));
        pst   = pst(q);
        y     = y(q) + Y(bin(q) + 1);

        % returned values
        %------------------------------------------------------------------
        if strcmp(xG.spec.Rplot,'fitted response and PSTH')
            G.PST  = PST;
            G.PSTH = PSTH;
            G.PCI  = PCI;
        end
        G.x   = x;
        G.CI  = CI;
        G.pst = pst;
        
        
    %-Parametric responses
    %======================================================================
    case 'Parametric responses'

        s     = xG.spec.Sess;
        u     = xG.spec.u;
        p     = xG.spec.p;
        
        % basis functions
        %------------------------------------------------------------------
        dt    = SPM.xBF.dt;
        bf    = SPM.xBF.bf;
        pst   = ([1:size(bf,1)] - 1)*dt;

        % orthogonalised expansion of parameteric variable
        %------------------------------------------------------------------
        P     = SPM.Sess(s).U(u).P(p).P;
        q     = [];
        for i = 0:SPM.Sess(s).U(u).P(p).h;
            q = [q P.^i];
        end
        q     = spm_orth(q);

        % parameter estimates for this effect
        %------------------------------------------------------------------
        B     = beta(SPM.Sess(s).col(SPM.Sess(s).Fc(u).i));

        % reconstruct trial-specific responses
        %------------------------------------------------------------------
        Y     = zeros(size(bf,1),size(q,1));
        uj    = SPM.Sess(s).U(u).P(p).i;
        for i = 1:size(P,1)
            U      = sparse(1,uj,q(i,:),1,size(SPM.Sess(s).U(u).u,2));
            X      = kron(U,bf);
            Y(:,i) = X*B;
        end
        [P,j] = sort(P);
        Y     = Y(:,j);

        % returned values
        %------------------------------------------------------------------
        G.pst = pst;
        G.P   = P;


    %-Volterra Kernels
    %======================================================================
    case 'Volterra Kernels'

        s     = xG.spec.Sess;
        u     = xG.spec.u;
        
        % Parameter estimates and basis functions
        %------------------------------------------------------------------
        dt    = SPM.xBF.dt;
        bf    = SPM.xBF.bf/dt;
        pst   = ([1:size(bf,1)] - 1)*dt;

        % second order kernel
        %------------------------------------------------------------------
        if u > length(SPM.Sess(s).U)

            % Parameter estimates and kernel
            %--------------------------------------------------------------
            B     = beta(SPM.Sess(s).col(SPM.Sess(s).Fc(u).i));
            i     = 1;
            Y     = 0;
            for p = 1:size(bf,2)
                for q = 1:size(bf,2)
                    Y = Y + B(i)*bf(:,p)*bf(:,q)';
                    i = i + 1;
                end
            end

        % first  order kernel
        %------------------------------------------------------------------
        else
            B = beta(SPM.Sess(s).col(SPM.Sess(s).Fc(u).i(1:size(bf,2))));
            Y = bf*B;
        end
        
        % returned values
        %------------------------------------------------------------------
        G.pst = pst;
        
end
