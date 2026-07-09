function D = bipolar_correct(D, flags)
%BIPOLAR_CORRECT Remove bipolar readout phase error (MEDI method) + noise map.
%
%   D = bipolar_correct(D, flags)
%
%   1. Removes the echo-dependent linear readout phase using the MEDI
%      iField_correction_new routine (bipolar / fly-back correction).
%   2. Estimates a noise standard-deviation map via Fit_ppm_complex.
%   3. Refines the tissue mask using that noise map.
%
%   flags fields (optional):
%       .verbose : logical, print progress (default true)
%       .nobone  : logical, remove bone from the refined mask (default false)
%
%   Condensed from BipolarCorrect (case 'MEDI') and the noise-STD block of the
%   original DogAnalysis.m. Only the MEDI correction is kept; the SEPIA and
%   Hernando bipolar variants were dropped as unused for this pipeline.

    if ~isfield(flags, 'verbose'), flags.verbose = true;  end
    if ~isfield(flags, 'nobone'),  flags.nobone  = false; end

    % ---- 1. Bipolar phase correction (MEDI) ----
    if ~D.Flags.CorrectedBipolarPhase
        if flags.verbose, fprintf('\nCorrecting for bipolar phase error...'); tic; end

        Mask = imdilate(D.Data.Mask, strel('disk', 4));
        D.Data.Image = iField_correction_new(D.Data.Image, D.VoxelSize, Mask, ...
                                              D.InPlanePhaseEncodingDirection);

        D.Flags.CorrectedBipolarPhase = true;
        D.BipolarCorrection = struct('Method', 'MEDI');
        if flags.verbose, tm = toc; fprintf('Done (%0.2f sec)', tm); end
    end

    % ---- 2. Noise standard-deviation map ----
    if ~isfield(D.Data, 'NoiseSTD')
        if flags.verbose, fprintf('\nCalculating noise STD...'); end
        opts.max_iter = 60;
        [~, D.Data.NoiseSTD, ~] = Fit_ppm_complex(D.Data.Image, opts);

        % ---- 3. Refine mask with noise map ----
        if flags.verbose, fprintf('\nRe-calculating binary mask with noise STD...'); end
        D.Data.Mask = GenMask(D.Data.Image, D.Data.WeightedMagnitude, flags.verbose, 3, ...
                              flags.nobone, D.Data.NoiseSTD);
    end
end
