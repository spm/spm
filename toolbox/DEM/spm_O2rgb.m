function I = spm_O2rgb(O,RGB)
% O to RGB image format
% FORMAT I = spm_O2rgb(O,RGB)
%
% Converts discrete (probabilistic) image into an RGB format for image
% display.
%
% This routine generates an image[s] using spatiotemporal basis functions
% weighted by expected singular variates; where the requisite expectations
% are evaluated over discrete values, given a multinomial distribution over
% those values. (If multiple images are supplied over the first is shown)
%
% see: spm_rgb2O
%__________________________________________________________________________


% if time series
%--------------------------------------------------------------------------
T   = size(O,2);
if T > 1
    I     = uint8(false(RGB.N(2),RGB.N(3),3,0));
    for t = 1:T
        I(:,:,:,end + (1:RGB.R)) = spm_O2rgb(O(:,t),RGB);
    end
    return
end

% Generate image
%--------------------------------------------------------------------------
I     = zeros(RGB.N);
for g = 1:numel(RGB.G)

    if isfield(RGB,'A')

        % recover singular vectors
        %------------------------------------------------------------------
        o     = find(ismember(RGB.O,g));
        Nm    = numel(o);
        u     = zeros(Nm,1);
        for m = 1:Nm
            if iscell(O)
                u(m) = RGB.A{o(m)}*O{o(m)};
            else
                u(m) = RGB.A{o(m)}(:,O(o(m)));
            end
        end
    else

        % or mutually exclusive images
        %------------------------------------------------------------------
        u     = double(O{g});
        
    end

    % place in image, unless zeros
    %----------------------------------------------------------------------
    if numel(u)
        I(RGB.G{g}) = I(RGB.G{g}) + double(RGB.V{g})*u;
    end

end

% permute colour to trailing dimension
%--------------------------------------------------------------------------
I = reshape(I,3,[],RGB.N(2),RGB.N(3));
I = uint8(permute(I,[3,4,1,2]));

return