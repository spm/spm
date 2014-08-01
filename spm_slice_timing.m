function spm_slice_timing(P, sliceorder, refslice, timing, prefix)
% Correct differences in slice acquisition times
% FORMAT spm_slice_timing(P, sliceorder, refslice, timing, prefix)
% P           - char array of image filenames
%               can also be a cell array of the above (multiple subjects).
% sliceorder  - slice acquisition order, a vector of integers, each
%               integer referring the slice number in the image file
%               (1=first), and the order of integers representing their
%               temporal acquisition order
%               OR vector containig the acquisition time for each slice
%               in milliseconds
% refslice    - slice for time 0
%               OR time in milliseconds for the reference slice
% timing      - additional information for sequence timing
%               timing(1) = time between slices
%                         = TA / (nslices - 1)
%               timing(2) = time between last slices and next volume
%                         = TR - TA
%               OR timing = [0 TR] when previous inputs are specified in
%               milliseconds
% prefix      - filename prefix for corrected image files, defaults to 'a'
%__________________________________________________________________________
%
%   Note: The sliceorder arg that specifies slice acquisition order is
%   a vector of N numbers, where N is the number of slices per volume.
%   Each number refers to the position of a slice within the image file.
%   The order of numbers within the vector is the temporal order in which
%   those slices were acquired.
%
%   To check the order of slices within an image file, use the SPM Display
%   option and move the crosshairs to a voxel co-ordinate of z=1.  This
%   corresponds to a point in the first slice of the volume.
%
%   The function corrects differences in slice acquisition times.
%   This routine is intended to correct for the staggered order of
%   slice acquisition that is used during echoplanar scanning. The
%   correction is necessary to make the data on each slice correspond
%   to the same point in time. Without correction, the data on one
%   slice will represent a point in time as far removed as 1/2 the TR
%   from an adjacent slice (in the case of an interleaved sequence).
%
%   This routine "shifts" a signal in time to provide an output
%   vector that represents the same (continuous) signal sampled
%   starting either later or earlier. This is accomplished by a simple
%   shift of the phase of the sines that make up the signal.
%
%   Recall that a Fourier transform allows for a representation of any
%   signal as the linear combination of sinusoids of different
%   frequencies and phases. Effectively, we will add a constant
%   to the phase of every frequency, shifting the data in time.
%
%   Shifter - This is the filter by which the signal will be convolved
%   to introduce the phase shift. It is constructed explicitly in
%   the Fourier domain. In the time domain, it may be described as
%   an impulse (delta function) that has been shifted in time the
%   amount described by TimeShift.
%
%   The correction works by lagging (shifting forward) the time-series
%   data on each slice using sinc-interpolation. This results in each
%   time series having the values that would have been obtained had
%   the slice been acquired at the same time as the reference slice.
%
%   To make this clear, consider a neural event (and ensuing hemodynamic
%   response) that occurs simultaneously on two adjacent slices. Values
%   from slice "A" are acquired starting at time zero, simultaneous to
%   the neural event, while values from slice "B" are acquired one
%   second later. Without corection, the "B" values will describe a
%   hemodynamic response that will appear to have began one second
%   EARLIER on the "B" slice than on slice "A". To correct for this,
%   the "B" values need to be shifted towards the Right, i.e., towards
%   the last value.
%
% Written by Darren Gitelman at Northwestern U., 1998
%
% Based (in large part) on ACQCORRECT.PRO from G. Aguirre and E. Zarahn
% at U. Penn.
%
% Modified by R. Henson, C. Buechel and J. Ashburner, FIL, to
% handle different reference slices and memory mapping.
%
% Modified by M. Erb, at U. Tuebingen, 1999, to ask for non-continuous
% slice timing and number of sessions.
%
% Modified by R. Henson for more general slice order and SPM2.
%
% Modified by A. Hoffmann, M. Woletz and C. Windischberger from Medical
% University of Vienna, Austria, to handle multi-band EPI sequences.
%__________________________________________________________________________
% Copyright (C) 1998-2014 Wellcome Trust Centre for Neuroimaging

% Darren Gitelman et al.
% $Id: spm_slice_timing.m 6130 2014-08-01 17:41:18Z guillaume $


SVNid = '$Rev: 6130 $';

%-Say hello
%--------------------------------------------------------------------------
SPMid = spm('FnBanner',mfilename,SVNid);

%-Parameters & Arguments
%==========================================================================
if nargin < 4, error('Not enough input arguments.'); end
if nargin < 5, prefix = 'a'; end

if ~iscell(P), P = {P}; end
nsubjects = numel(P);

% Acquisition order: 1=first slice in image
% Reference slice: 1=first slice in image, in Analyze format, slice 1 = bottom
% TR: Interscan interval (TR) {secs}
% TA: Acquisition Time (TA) {secs} [Def: TR-TR/nslices], TA <= TR
% timing(2) = TR - TA, time between last slices and next volume
% timing(1) = TA / (nslices -1), time between slices

Vin     = spm_vol(P{1}(1,:));
nslices = Vin(1).dim(3);

TR      = (nslices-1)*timing(1)+timing(2);
fprintf('%-40s: %30s\n','Number of slices is...',num2str(nslices))      %-#
fprintf('%-40s: %30s\n','Time to Repeat (TR) is...',num2str(TR))        %-#

if ~isequal(1:nslices,sort(sliceorder))
    if ~all(sliceorder >= 0 & sliceorder <= TR*1000)
        error('Input is neither slice indices nor slice times.');
    end
    unit = 'slice times (ms)';
else
    if ~ismember(refslice,sliceorder)
        error('Reference slice should contain a slice index.');
    end
    unit = 'slice indices';
end
fprintf('%-40s: %30s\n','Parameters are specified as...',unit)          %-#

if nslices ~= numel(sliceorder)
    error('Mismatch between number of slices and length of ''sliceorder'' vector.');
end

%-Slice timing correction
%==========================================================================
for subj = 1:nsubjects
    Vin   = spm_vol(P{subj});
    nimgo = numel(Vin);
    nimg  = 2^(floor(log2(nimgo))+1);
    if Vin(1).dim(3) ~= nslices
        error('Number of slices differ: %d vs %d.', nslices, Vin(1).dim(3));
    end
        
    % Create new header files
    Vout  = Vin;
    for k=1:nimgo
        Vout(k).fname  = spm_file(Vin(k).fname, 'prefix', prefix);
        if isfield(Vout(k),'descrip')
            desc = [Vout(k).descrip ' '];
        else
            desc = '';
        end
        Vout(k).descrip = [desc 'acq-fix ref-slice ' num2str(refslice)];
    end
    Vout = spm_create_vol(Vout);
    
    % Set up [time x voxels] matrix for holding image info
    slices = zeros([Vout(1).dim(1:2) nimgo]);
    stack  = zeros([nimg Vout(1).dim(1)]);
    
    task = sprintf('Correcting acquisition delay: session %d', subj);
    spm_progress_bar('Init',nslices,task,'planes complete');
    
    % Compute shifting amount from reference slice and slice order
    if isequal(unit,'slice times (ms)')
        % Compute time difference between the acquisition time of the
        % reference slice and the current slice by using slice times
        % supplied in sliceorder vector
        shiftamount = (sliceorder - refslice) / (1000 * TR);
    else
        rslice      = find(sliceorder==refslice);
        [Y, I]      = sort(sliceorder);
        shiftamount = (I - rslice) * timing(1) / TR;
    end
    
    % For loop to perform correction slice by slice
    for k = 1:nslices
        
        % Read in slice data
        B  = spm_matrix([0 0 k]);
        for m=1:nimgo
            slices(:,:,m) = spm_slice_vol(Vin(m),B,Vin(1).dim(1:2),1);
        end
        
        % Set up shifting variables
        len     = size(stack,1);
        phi     = zeros(1,len);
        
        % Check if signal is odd or even -- impacts how Phi is reflected
        %  across the Nyquist frequency. Opposite to use in pvwave.
        OffSet  = 0;
        if rem(len,2) ~= 0, OffSet = 1; end
        
        % Phi represents a range of phases up to the Nyquist frequency
        % Shifted phi 1 to right.
        for f = 1:len/2
            phi(f+1) = -1*shiftamount(k)*2*pi/(len/f);
        end
        
        % Mirror phi about the center
        % 1 is added on both sides to reflect Matlab's 1 based indices
        % Offset is opposite to program in pvwave again because indices are 1 based
        phi(len/2+1+1-OffSet:len) = -fliplr(phi(1+1:len/2+OffSet));
        
        % Transform phi to the frequency domain and take the complex transpose
        shifter = [cos(phi) + sin(phi)*sqrt(-1)].';
        shifter = shifter(:,ones(size(stack,2),1)); % Tony's trick
        
        % Loop over columns
        for i=1:Vout(1).dim(2)
            
            % Extract columns from slices
            stack(1:nimgo,:) = reshape(slices(:,i,:),[Vout(1).dim(1) nimgo])';
            
            % Fill in continous function to avoid edge effects
            for g=1:size(stack,2)
                stack(nimgo+1:end,g) = linspace(stack(nimgo,g),...
                    stack(1,g),nimg-nimgo)';
            end
            
            % Shift the columns
            stack = real(ifft(fft(stack,[],1).*shifter,[],1));
            
            % Re-insert shifted columns
            slices(:,i,:) = reshape(stack(1:nimgo,:)',[Vout(1).dim(1) 1 nimgo]);
        end
        
        % Write out the slice for all volumes
        for p = 1:nimgo
            Vout(p) = spm_write_plane(Vout(p),slices(:,:,p),k);
        end
        spm_progress_bar('Set',k);
    end
    spm_progress_bar('Clear');
end

fprintf('%-40s: %30s\n','Completed',spm('time'))                        %-#
