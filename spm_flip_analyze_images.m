% Do Analyze format images need to be left-right flipped?
%-----------------------------------------------------------------------
function flip = spm_flip_analyze_images


%========================== START DELETE ===============================
% Issue an error if not already configured..
%-----------------------------------------------------------------------
if exist('isFIL') ~= 2,
	error('You need to configure this function for your particular setup.\n');
end;
%============================ END DELETE ===============================


% Uncomment the relevant lines below...
%-----------------------------------------------------------------------
flip = 1; % Default Analyze orientation
% flip = 0; % Left-right flipped Analyze orientation
