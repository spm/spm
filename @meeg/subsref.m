function varargout=subsref(this,subs)
% SUBSREF Subscripted reference
% An overloaded function...
% _________________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak, Stefan Kiebel
% $Id: subsref.m 1243 2008-03-25 23:02:44Z stefan $

if isempty(subs)
    return;
end;

if this.Nsamples == 0
    error('Attempt to reference a field of an empty meeg object.');
end;

switch subs(1).type
    case '()'
        if numel(subs)~= 1, error('Expression too complicated');end;
        varargout = {subsref(this.data.y, subs)};
    case '{}'
    case '.'
        if isfield(this.other, subs(1).subs)
            field = getfield(this.other, subs(1).subs);
            if numel(subs)==1
                varargout = {field};
            else
                varargout ={subsref(field, subs(2:end))};
            end
        else
            switch(subs(1).subs)
                case 'nchannels'
                    varargout = {length(this.channels)};
                case 'nconditions'
                    varargout = {size(unique(conditions(this), 'rows'),1)};
                case 'nsamples'
                    varargout = {this.Nsamples};
                case 'dtype';
                    % returns datatype of embedded file_array object
                    varargout = {this.data.y.dtype};
                case 'fsample'
                    varargout = {this.Fsample};
                case 'meegchannels'
                    type = cat(1, this.channels(:).type);
                    varargout = {unique([find(strmatch('EEG', type)) find(strmatch('MEG', type))])};
                case 'ntrials'
                    varargout = {length(this.trials)};

                otherwise
                    if ismethod(this, subs(1).subs)
                        switch numel(subs)
                            case 1
                                varargout = {feval(subs(1).subs, this)};
                            case 2
                                varargout = {feval(subs(1).subs, this, subs(2).subs{:})};
                            otherwise
                                error('Expression too complicated');
                        end
                    else
                        error('Reference to non-existent or private meeg method or field.');
                    end
            end
        end
    otherwise
        error('Unfamiliar referencing type');
end


