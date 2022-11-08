function spm_opm_toot(snd)
% Audibly notify the end of a script, or the working day perhaps?
% FORMAT spm_opm_toot(snd)
% snd   - name of file containing a sound [Default: 'train']
%         or one of ['chirp','gong','handel','laughter','splat','train']
%__________________________________________________________________________

% George O'Neill
% Copyright (C) 2022 Wellcome Centre for Human Neuroimaging


if ~nargin, snd = 'train'; end
try
    snd = load(snd);
    sound(snd.y,snd.Fs);
end
