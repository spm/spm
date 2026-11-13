function spm_imshow(I)
% rendering of an RGB image sequence
% FORMAT spm_imshow(I)
%
% I(pixels x pixels x rgb x time)
%__________________________________________________________________________

% image
%--------------------------------------------------------------------------
for t = 1:size(I,4)
   imshow(uint8(I(:,:,:,t)))
   drawnow
end

return