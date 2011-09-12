function spm_slice_timing(P, sliceorder, refslice, timing, prefix)
% Correct differences in slice acquisition times
% FORMAT spm_slice_timing(P, sliceorder, refslice, timing, prefix)
% P           - nimages x ? Matrix with filenames
%               can also be a cell array of the above (multiple subj).
% sliceorder  - slice acquisition order, a vector of integers, each
%               integer referring the slice number in the image file
%               (1=first), and the order of integers representing their
%               temporal acquisition order
% refslice    - slice for time 0
% timing      - additional information for sequence timing
%               timing(1) = time between slices
%               timing(2) = time between last slices and next volume
% prefix      - filename prefix for corrected image files, defaults to 'a'
%
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
%   This correction assumes that the data are band-limited (i.e. there
%   is no meaningful information present in the data at a frequency
%   higher than that of the Nyquist). This assumption is support by
%   the study of Josephs et al (1997, NeuroImage) that obtained
%   event-related data at an effective TR of 166 msecs. No physio-
%   logical signal change was present at frequencies higher than our
%   typical Nyquist (0.25 HZ).
%
% Written by Darren Gitelman at Northwestern U., 1998
%
% Based (in large part) on ACQCORRECT.PRO from Geoff Aguirre and
% Eric Zarahn at U. Penn.
%
% v1.0  07/04/98    DRG
% v1.1  07/09/98    DRG fixed code to reflect 1-based indices
%               of matlab vs. 0-based of pvwave
%
% Modified by R Henson, C Buechel and J Ashburner, FIL, to
% handle different reference slices and memory mapping.
%
% Modified by M Erb, at U. Tuebingen, 1999, to ask for non-continuous
% slice timing and number of sessions.
%
% Modified by R Henson for more general slice order and SPM2
%__________________________________________________________________________
% Copyright (C) 1999-2011 Wellcome Trust Centre for Neuroimaging

% Darren Gitelman et al.
% $Id: spm_slice_timing.m 4479 2011-09-12 11:28:04Z guillaume $


SVNid = '$Rev: 4479 $';

%-Say hello
%--------------------------------------------------------------------------
SPMid = spm('FnBanner',mfilename,SVNid);

%-Parameters & Arguments
%==========================================================================
if nargin < 4
    error('Not enough input arguments.');
end

if iscell(P)
    nsubjects = length(P);
else
    nsubjects = 1;
    P = {P};
end

% Acquisition order: 1=first slice in image
% Reference slice: 1=first slice in image, in Analyze format, slice 1 = bottom
% TR: Interscan interval (TR) {secs}
% TA: Acquisition Time (TA) {secs} [Def: TR-TR/nslices], TA <= TR
% timing(2) = TR - TA, time between slices
% timing(1) = TA / (nslices -1), time between last slices and next volume

Vin     = spm_vol(P{1}(1,:));
nslices = Vin(1).dim(3);

TR  = (nslices-1)*timing(1)+timing(2);
fprintf('%-40s: %30s\n','Number of slices is...',num2str(nslices))      %-#
fprintf('%-40s: %30s\n','Time to Repeat (TR) is...',num2str(TR))        %-#
factor = timing(1)/TR;

if nargin < 5, prefix = 'a'; end

%-Slice timing correction
%==========================================================================
for subj = 1:nsubjects
    PP        = P{subj};
    Vin       = spm_vol(PP);
    nimgo     = numel(Vin);
    nimg      = 2^(floor(log2(nimgo))+1);
    if Vin(1).dim(3) ~= nslices
        error('Number of slices differ! %d %\n', nimg);
    end
    
    % create new header files
    Vout    = Vin;
    for k=1:nimgo
        Vout(k).fname  = spm_file(Vin(k).fname, 'prefix', prefix);
        if isfield(Vout(k),'descrip')
            desc = [Vout(k).descrip ' '];
        else
            desc = '';
        end
        Vout(k).descrip = [desc 'acq-fix ref-slice ' int2str(refslice)];
    end
    Vout = spm_create_vol(Vout);
    
    % Set up large matrix for holding image info
    % Organization is time by voxels
    slices = zeros([Vout(1).dim(1:2) nimgo]);
    stack  = zeros([nimg Vout(1).dim(1)]);
    
    task = sprintf('Correcting acquisition delay: session %d', subj);
    spm_progress_bar('Init',nslices,task,'planes complete');
    
    % For loop to read data slice by slice do correction and write out
    % In analzye format, the first slice in is the first one in the volume.
    
    rslice = find(sliceorder==refslice);
    for k = 1:nslices
        
        % Set up time acquired within slice order
        shiftamount  = (find(sliceorder==k) - rslice) * factor;
        
        % Read in slice data
        B  = spm_matrix([0 0 k]);
        for m=1:nimgo
            slices(:,:,m) = spm_slice_vol(Vin(m),B,Vin(1).dim(1:2),1);
        end
        
        % set up shifting variables
        len     = size(stack,1);
        phi     = zeros(1,len);
        
        % Check if signal is odd or even -- impacts how Phi is reflected
        %  across the Nyquist frequency. Opposite to use in pvwave.
        OffSet  = 0;
        if rem(len,2) ~= 0, OffSet = 1; end
        
        % Phi represents a range of phases up to the Nyquist frequency
        % Shifted phi 1 to right.
        for f = 1:len/2
            phi(f+1) = -1*shiftamount*2*pi/(len/f);
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
            
            % fill in continous function to avoid edge effects
            for g=1:size(stack,2)
                stack(nimgo+1:end,g) = linspace(stack(nimgo,g),...
                    stack(1,g),nimg-nimgo)';
            end
            
            % shift the columns
            stack = real(ifft(fft(stack,[],1).*shifter,[],1));
            
            % Re-insert shifted columns
            slices(:,i,:) = reshape(stack(1:nimgo,:)',[Vout(1).dim(1) 1 nimgo]);
        end
        
        % write out the slice for all volumes
        for p = 1:nimgo
            Vout(p) = spm_write_plane(Vout(p),slices(:,:,p),k);
        end
        spm_progress_bar('Set',k);
    end
    spm_progress_bar('Clear');
end

fprintf('%-40s: %30s\n','Completed',spm('time'))                        %-#
