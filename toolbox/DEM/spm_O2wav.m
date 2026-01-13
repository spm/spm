function wav = spm_O2wav(O,WAV)
% O to WAV image format
% FORMAT wav = spm_O2wav(O,WAV)
%
% Converts discrete (probabilistic) image into an WAV format for sonogram
% [dis]play.
%
% This routine generates an image using spatiotemporal basis functions
% weighted by expected singular variance; where the requisite expectations
% are evaluated over discrete values, given a multinomial distribution over
% those values. This routine is specialised for one-dimensional timeseries;
% e.g., time frequency or sonogram representation of audio files.
%
% see also: spm_wav2O and spm_O2rgb
%__________________________________________________________________________


% Generate sound file
%--------------------------------------------------------------------------
T     = size(O,2);
wav   = cell(1,T);
for t = 1:T
    I     = zeros(WAV.N);
    for g = 1:numel(WAV.G)

        if isfield(WAV,'A')

            % recover singular vectors
            %--------------------------------------------------------------
            o     = find(ismember(WAV.O,g));
            Nm    = numel(o);
            u     = zeros(Nm,1);
            for m = 1:Nm
                if iscell(O)
                    u(m) = WAV.A{o(m)}*O{o(m),t};
                else
                    u(m) = WAV.A{o(m)}(:,O(o(m),t));
                end
            end
        else

            % or mutually exclusive images
            %--------------------------------------------------------------
            u     = double(O{g});

        end

        % place in image, unless zeros
        %------------------------------------------------------------------
        if numel(u)
            I(WAV.G{g}) = I(WAV.G{g}) + double(WAV.V{g})*u;
        end

    end

    % store this segement
    %----------------------------------------------------------------------
    wav{1,t} = I';

end

% concatentate segements
%--------------------------------------------------------------------------
wav = spm_cat(wav);

return