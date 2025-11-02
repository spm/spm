function spm_bsd_results(BSD, fignum)
% spm_bsd_results - Processes BSD results and optionally displays figures.
%
% Syntax:
%   spm_bsd_results(BSD, fignum)
%
% Inputs:
%   BSD    - Structure containing BSD analysis results.
%   fignum - (Optional) Figure number for result visualization.
%
% Description:
%   This function processes the results contained in the BSD structure,
%   performing necessary computations and optionally displaying results
%   in the specified figure window.
%_______________________________________________________________________
% Copyright (C) 2024-2025 Wellcome Trust Centre for Neuroimaging

% Johan Medrano
    
if nargin < 2
    fignum = 1; 
end

if iscell(BSD)
    for i = 1:numel(BSD)
        spm_bsd_results(BSD{i}, i);
    end
    return 
elseif (isstruct(BSD) && isfield(BSD, 'models'))
    for i = 1:numel(BSD.models)
        spm_bsd_results(BSD.models(i), i);
    end
    return 
else
    bsd = BSD;
end


% create spm figure
%--------------------------------------------------------------------------
spm_figure('GetWin', sprintf('BSD - Results (%d)', fignum));
clf;

% apply exponential transform if in log
%--------------------------------------------------------------------------
if bsd.M.fitlog
    trf  = @exp; 
    trfi = @log; 
else 
    trf = @(x) real(x); 
    trfi = trf; 
end

%% plot credible intervals
%==========================================================================

% use upper half
%--------------------------------------------------------------------------
subplot(2,1,1);
cla

% prepare precision components
%--------------------------------------------------------------------------
Q = bsd.xY.Q;
Q = bsd.Ce * (1./(diag(Q)+1e-4)); 

% compute 95% credible intervals
%--------------------------------------------------------------------------
m = bsd.Hs{1}; 
s = sqrt(Q);  
l = m - 1.96 * s;
u = m + 1.96 * s;

% check if powerline was provided... 
%--------------------------------------------------------------------------
if isempty(bsd.M.powerline) 
    % ... if absent, plot full intervals
    %----------------------------------------------------------------------
    fqs = bsd.xY.Hz; 
    fill([fqs, flip(fqs)]', trf([u; flip(l)]), [.6,.6,.6], ...
         'FaceAlpha', 0.3, 'EdgeColor', 'white')
    hold on; 
else 
    % ... otherwise, separate the powerline interval
    %======================================================================

    % get powerline 
    %----------------------------------------------------------------------
    pwl = bsd.M.powerline; 
    
    % plot left interval
    %----------------------------------------------------------------------
    fsel = bsd.xY.Hz <= pwl(1); 
    fqs  = bsd.xY.Hz(fsel); 
    uL   = u(fsel); 
    lL   = l(fsel);
    fill([fqs, flip(fqs)]', trf([uL; flip(lL)]), [.6,.6,.6], ...
        'FaceAlpha', 0.1, 'EdgeAlpha', 0)
    hold on; 

    % plot powerline
    %----------------------------------------------------------------------
    fsel = bsd.xY.Hz >= pwl(1) & bsd.xY.Hz <= pwl(2); 
    fqs  = bsd.xY.Hz(fsel); 
    uP   = u(fsel); 
    lP   = l(fsel);
    fill([fqs, flip(fqs)]', trf([uP; flip(lP)]), 'r', ...
        'FaceAlpha', 0.1, 'EdgeAlpha', 0)

    % plot right interval
    %----------------------------------------------------------------------
    fsel = bsd.xY.Hz >= pwl(2); 
    fqs  = bsd.xY.Hz(fsel); 
    uR   = u(fsel); 
    lR   = l(fsel);
    fill([fqs, flip(fqs)]', trf([uR; flip(lR)]), [.6,.6,.6], ...
        'FaceAlpha', 0.1, 'EdgeAlpha', 0)
end

% plot periodic components
%==========================================================================

% get modes parameters
%--------------------------------------------------------------------------
spec = spm_bsd_param2spec(bsd.Ep, bsd.M);

% plot 1/f
%--------------------------------------------------------------------------
plot(bsd.xY.Hz, trf(bsd.Hsens{1}.noise), 'k')

% plot modes
%--------------------------------------------------------------------------
colors = colororder; 
for i = 1:numel(bsd.fqs)
    % extract frequency band
    %----------------------------------------------------------------------
    fsel = bsd.xY.Hz >= bsd.fqs{i}(1) & bsd.xY.Hz <= bsd.fqs{i}(2); 
    fqs  = bsd.xY.Hz(fsel);
    ampl = bsd.Hsens{1}.modes{i}(fsel); 
    grnd = bsd.Hsens{1}.noise(fsel); 
    
    % show mode
    %----------------------------------------------------------------------
    c = colors(mod(i, length(colors)), :); % mode color
    fill([fqs flip(fqs)]', trf([grnd; flip(grnd+ampl)]), c, ...
        'FaceAlpha', 0.3, 'EdgeAlpha', 0);
    plot(fqs, trf(grnd+ampl), '-', 'Color', c, 'LineWidth', 2)
    
    % annotate mode
    %----------------------------------------------------------------------
    xtext = [spec.freq(i) spec.freq(i) + 1]; 
    ytext = interp1(fqs, grnd, spec.freq(i), 'cubic') + spec.ampl(i) + [0 1];
    str   = sprintf('%.1fHz', spec.freq(i)); 
    text(xtext(2), trf(ytext(2)), str, 'Color',  0.2*c);
    plot(xtext, trf(ytext), '-', 'LineWidth', 1.1, 'Color', 0.2*c);
end

% plot observed data
%--------------------------------------------------------------------------
plot(bsd.xY.Hz,  trf(bsd.xY.y{1}), 'k:', 'LineWidth', 1.5)

% finalize plot 
%==========================================================================

% setup axes ranges
%--------------------------------------------------------------------------
ymax = max([bsd.Hsens{1}.modes{:}] + bsd.Hsens{1}.noise, [], 'all');
ymax = max(ymax, max(bsd.xY.y{1}));
ymin = min([bsd.Hsens{1}.modes{:}] + bsd.Hsens{1}.noise, [], 'all');
ymin = min(ymin, min(bsd.xY.y{1}));
axis tight
ylim([trf(ymin) trf(ymin + 1.2 * (ymax-ymin))])

% annotate powerline
%--------------------------------------------------------------------------
pwl = bsd.M.powerline; 
text((pwl(1)+pwl(2))/2+0.1, trf(ymin + 1.1 * (ymax-ymin)), 'powerline', ...
    'Rotation', -90, 'Color', [0.6 0 0], 'VerticalAlignment', 'middle')

% setup axes titles
%--------------------------------------------------------------------------
xlabel('Frequency (Hz)')
if bsd.M.fitlog
    ylabel('Log amplitude (a.u.)');
    yscale('log');
else 
    ylabel('Amplitude (a.u)');
end

% setup plot title
%--------------------------------------------------------------------------
title('Bayesian spectral decomposition')

%% plot components
%==========================================================================

% helper to compute variance using delta transformation
%--------------------------------------------------------------------------
specvar = @(E, C, M, spec) spm_unvec( ...
    spm_cat(spm_diff(@spm_bsd_param2spec, E, M, 1))  ...
  * spm_vec(C) ...
  , spec);

% priors in specification space
%--------------------------------------------------------------------------
sE = spm_bsd_param2spec(bsd.M.pE, bsd.M);
sC = specvar(bsd.M.pE, bsd.M.pC, bsd.M, sE); 

% posteriors in specification space
%--------------------------------------------------------------------------
Es = spm_bsd_param2spec(bsd.Ep, bsd.M);
Cs = specvar(bsd.Ep, diag(bsd.Cp), bsd.M, Es); 

% handle independent fit of aperiodic component
%--------------------------------------------------------------------------
if bsd.M.separatenull
    % prior in specification space (from null model)
    %----------------------------------------------------------------------
    rE = spm_bsd_param2spec(bsd.null.M.pE, bsd.null.M);
    rC = specvar(bsd.null.M.pE, bsd.null.M.pC, bsd.null.M, rE); 

    sE.noise = rE.noise; 
    sC.noise = rC.noise; 

    % posterior in specification space (from null model)
    %----------------------------------------------------------------------
    rE = spm_bsd_param2spec(bsd.null.Ep, bsd.null.M);
    rC = specvar(bsd.null.Ep, diag(bsd.null.Cp), bsd.null.M, rE); 

    Es.noise = rE.noise; 
    Cs.noise = rC.noise; 
end

% prepare labels and titles for plotting 
%--------------------------------------------------------------------------
fields = {'freq', 'fwhm', 'ampl', '', '', ''};
titles = {'Frequency', 'FWHM', 'Amplitude', ...
    'Intercept', 'Exponent', 'Knee frequency'}; 
ylabels = {'Frequency (Hz)', 'Frequency (Hz)', '', '', ...
    'Frequency (Hz)', 'Frequency (Hz)'}; 
if bsd.M.fitlog
    ylabels{3} = 'Log amplitude (a.u.)'; 
    ylabels{4} = 'Log amplitude (a.u.)'; 
else
    ylabels{3} = 'Amplitude (a.u)'; 
    ylabels{4} = 'Amplitude (a.u)'; 
end

% Iterate over specification parameters
%--------------------------------------------------------------------------
for i = 1:6
    % setup plot
    %----------------------------------------------------------------------
    subplot(4, 3, 6+i)
    title(titles{i}); 
    hold on; 

    % setup scale
    %----------------------------------------------------------------------
    if (i == 3 || i == 4) && bsd.M.fitlog
        yscale log
        scale = 'log'; 
    end
    
    % get fieldname 
    %----------------------------------------------------------------------
    f = fields{i}; 

    % compute prior and posterior expectation and covariance...
    %----------------------------------------------------------------------
    if i < 3 
        % ... for frequency and fwhm
        %------------------------------------------------------------------
        pE = [sE(:).(f)]; 
        pC = 1.96 * sqrt([sC(:).(f)]); 
        Ep = [Es(:).(f)];
        Cp = 1.96 * sqrt([Cs(:).(f)]); 

    elseif i == 3
        % ... for mode amplitude
        %------------------------------------------------------------------
        pE = trfi([sE(:).(f)]); 
        pC = 1.96 * sqrt(bsd.M.pC.a); 
        Ep = trfi([Es(:).(f)]);
        Cp = spm_unvec(1.96 * sqrt(diag(bsd.Cp)), bsd.Ep);
        Cp = Cp.a;

    elseif i == 4
        % ... for aperiodic intercept
        %------------------------------------------------------------------
        pE = trfi([sE(:).noise(1)]); 
        pC = 1.96 * sqrt(bsd.null.M.pC.b(1)); 
        Ep = trfi([Es(:).noise(1)]);
        Cp = spm_unvec(1.96 * sqrt(diag(bsd.null.Cp)), bsd.null.Ep);
        Cp = Cp.b(1); 

    elseif i == 5 
        % ... for aperiodic exponent 
        %------------------------------------------------------------------
        pE = [sE.noise(2)]; 
        pC = 1.96 * sqrt([sC.noise(2)]); 
        Ep = [Es.noise(2)];
        Cp = 1.96 * sqrt([Cs.noise(2)]);

    elseif i == 6 
        % ... for knee frequency (using delta method)
        %==================================================================

        % helper conversion function
        %------------------------------------------------------------------
        kEf = @(k, e) k.^(1/e); 
        kCf = @(k, kC, e, eC) ...
            spm_diff(kEf, k, e, 1) * kC + spm_diff(kEf, k, e, 2) * eC; 

        % prior moments (parameter space)
        %------------------------------------------------------------------
        kE = sE.noise(3); % prior expectation - knee parameter
        kC = sC.noise(3); % prior variance - knee parameter
        eE = sE.noise(2); % prior expectation - exponent parameter
        eC = sC.noise(2); % prior variance - exponent parameter
        
        % prior moments for knee frequency
        %------------------------------------------------------------------
        pE = kEf(kE, eE);         % prior expectation
        pC = kCf(kE, kC, eE, eC); % prior variance

        % posterior moments (parameter space)
        %------------------------------------------------------------------
        Ek = Es.noise(3); % posterior expectation - knee parameter
        Ck = Cs.noise(3); % posterior variance - knee parameter
        Ee = Es.noise(2); % posterior expectation - exponent parameter
        Ce = Cs.noise(2); % posterior variance - exponent parameter
        
        % posterior moments for knee frequency
        %------------------------------------------------------------------
        pE = kEf(Ek, Ee);         % prior expectation
        pC = kCf(Ek, Ck, Ee, Ce); % prior variance
    end

    % construct plot
    %----------------------------------------------------------------------
    dx = 0.1;
    for j = 1:length(pE)
        % setup indices
        %------------------------------------------------------------------
        x = [j-dx j+dx]; 
        yE = [pE(j) Ep(j)]; 
        yC = [pC(j) Cp(j)]; 
        
        % compute intervals
        %------------------------------------------------------------------
        if (i == 3 || i == 4) && bsd.M.fitlog % handle log fit
            yP = (exp(+yC) - 1).*exp(yE);
            yM = (exp(-yC) - 1).*exp(yE);
            yL = exp(yE - yC); 
            yU = exp(yE + yC); 
            yE = exp(yE);
        else
            yP = yC; 
            yM = -yC; 
            yL = yE + yM; 
            yU = yE + yP; 
        end

        % plot credible interval block
        %------------------------------------------------------------------
        c = colors(j, :);
        fill([x flip(x)], [yL, flip(yU)], c, ... 
            'FaceAlpha', 0.2, 'EdgeColor', c)

        % link prior and posterior expectations
        %------------------------------------------------------------------
        plot(x, yE, 'w', 'LineWidth', 2)

        % plot prior error bars 
        %------------------------------------------------------------------
        c = 0.6*[1 1 1];
        errorbar(j-dx, yE(1), yM(1), yP(1), 'Marker', 'o', 'Color', c, ...
            'MarkerFaceColor', 'w', 'MarkerEdgeColor', c, ...
            'MarkerSize', 10, 'LineWidth', 1.5, 'LineStyle', '-'); 

        % plot posterior error bars
        %------------------------------------------------------------------
        c = colors(j, :);
        errorbar(j+dx, yE(2), yM(2), yP(2), 'Marker', 'o', 'Color', c, ...
            'MarkerFaceColor', 'w', 'MarkerEdgeColor', c, ...
            'MarkerSize', 10, 'LineWidth', 1.5, 'LineStyle', '-');
    end
    
    % setup axes and legend
    %----------------------------------------------------------------------
    ylabel(ylabels{i});
    xlim([0 length(pE)+1]);

    lims = ylim; 
    lims(2) = lims(1) + 1.2 * (lims(2) - lims(1)); 
    ylim(lims);
    legend({'', '', 'prior', 'posterior'}, 'Box', 'off', ...
        'Orientation', 'horizontal', 'Location', 'north')

end 