function res = type(this, value)
% Method for and getting/setting EEG file type
% FORMAT res = type(this, value)
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


if nargin == 1
    res = this.type;
else
    switch value
        case 'continuous'
            if ntrials(this)>1
                error('Continuous file can only have one trial');
            end
        case 'single'
        case 'evoked' % Add additional checks here
        case 'grandmean'
        otherwise
            error('Unrecognized type');
    end

    this.type = value;
    res = this;
end
